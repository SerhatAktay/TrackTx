import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.MessageDigest

/**
 * Stage store: content-addressed result directories for every pipeline stage.
 *
 * Each process sets `storeDir { Stamp.dir(...) }`. The directory name embeds a
 * hash (the "stamp") of everything that can change the stage's output:
 *   - the module source and every bin/ script it calls (found by scanning the text)
 *   - the params the module references (minus the ones that never change results)
 *   - the container tag and lib/*.groovy
 *   - the stage inputs: val values as-is; files by store-relative path when they
 *     come from another stage (which already carries that stage's stamp, so
 *     stamps chain), else by content (small) or name+size (large)
 *
 * If all declared outputs already exist in that directory Nextflow skips the task
 * and re-emits them, with or without the work/ dir. Change a module, a script or a
 * param and only the affected stage's stamp changes; every stage downstream of it
 * changes too because its input paths change. A human-readable record of each
 * stamp is written to <output_dir>/.store/_manifests/<module>/<label>-<stamp>.txt.
 *
 * NOTE: storeDir closures see params from -params-file / nextflow.config, but NOT
 * undeclared command-line `--key value` overrides (Nextflow limitation).
 */
class Stamp {

  // Params that never change a stage's scientific output. Everything else a module
  // references is hashed, so a missing entry here only costs a needless rerun.
  static final Set<String> IGNORE = [
    'output_dir', 'publish_mode', 'publish_trimmed_fastq', 'publish_alignments',
    'publish_sra_fastq', 'publish_references_gtf', 'verbose', 'debug', 'help',
    'samplesheet', 'experiment_name', 'assets_dir', 'genome_cache', 'sra_cache_dir',
    'sra_tmp', 'sra_max_size', 'sra_download_source', 'conda_pol', 'conda_norm',
    'conda_fgr', 'conda_divergent', 'conda_annotation', 'conda_sra', 'get'
  ] as Set

  static final int INLINE_MAX = 8 * 1024 * 1024   // hash content up to this size

  private static final Map<String, String> TEXT = [:]

  static String dir(String module, Map params, Path projectDir, List inputs, String label) {
    def root  = new File(params.output_dir.toString()).absoluteFile.toPath().normalize()
    def store = root.resolve('.store')
    def parts = signature(module, params, projectDir, inputs, root)
    def stamp = sha(parts.join('\n')).take(16)
    def name  = "${label}-${stamp}".toString().replaceAll(/[^A-Za-z0-9_.\-]/, '_')
    writeManifest(store, module, name, parts)
    return store.resolve(module).resolve(name).toString()
  }

  static List<String> signature(String module, Map params, Path projectDir, List inputs, Path root) {
    def parts = []
    def modFile = new File(projectDir.toString(), "modules/${module}.nf")
    def modText = read(modFile)
    parts << "module=${module} code=${sha(modText)}".toString()
    deps(projectDir, modText).each { n -> parts << "bin/${n}=${sha(read(new File(projectDir.toString(), "bin/${n}")))}".toString() }
    parts << "lib=${sha(libText(projectDir))}".toString()
    parts << "container=${containerTag(projectDir)}".toString()
    // containsKey: a name found only in module comments (e.g. "params.yaml") is not a
    // real param, and reading it with [] would make Nextflow warn about an undefined one.
    paramKeys(modText).each { k -> parts << "param.${k}=${canon(params.containsKey(k) ? params[k] : null, root)}".toString() }
    inputs.eachWithIndex { v, i -> parts << "in${i}=${canon(v, root)}".toString() }
    return parts
  }

  // ── source scanning ──────────────────────────────────────────────────────

  private static synchronized String read(File f) {
    TEXT.computeIfAbsent(f.path) { p -> f.exists() ? f.text : '' }
  }

  private static String libText(Path projectDir) {
    def d = new File(projectDir.toString(), 'lib')
    (d.listFiles() ?: [] as File[]).findAll { it.name.endsWith('.groovy') }.sort { it.name }
      .collect { read(it) }.join('\n')
  }

  private static String containerTag(Path projectDir) {
    def m = (read(new File(projectDir.toString(), 'nextflow.config')) =~ /ghcr\.io\/[\w.\-\/]+:([\w.\-]+)/)
    m.find() ? m.group(1) : 'unknown'
  }

  // bin/ scripts named in the module text, plus anything they call in turn.
  private static synchronized List<String> deps(Path projectDir, String modText) {
    def binDir = new File(projectDir.toString(), 'bin')
    def names  = (binDir.listFiles() ?: [] as File[]).findAll { it.isFile() }*.name as Set
    def seen = new TreeSet<String>()
    def todo = [modText]
    while (!todo.isEmpty()) {
      def text = todo.remove(0)
      names.each { n ->
        if (!seen.contains(n) && text.contains(n)) {
          seen << n
          todo << read(new File(binDir, n))
        }
      }
    }
    return seen.toList()
  }

  private static List<String> paramKeys(String modText) {
    def keys = new TreeSet<String>()
    def m = (modText =~ /params\.(?:get\(\s*['"]([A-Za-z_]\w*)['"]|([A-Za-z_]\w*))/)
    while (m.find()) {
      def k = m.group(1) ?: m.group(2)
      if (!IGNORE.contains(k)) keys << k
    }
    return keys.toList()
  }

  // ── canonical values ─────────────────────────────────────────────────────

  static String canon(Object v, Path root) {
    if (v == null) return 'null'
    if (v instanceof Path) return fingerprint(v, root)
    if (v instanceof File) return fingerprint(v.toPath(), root)
    if (v instanceof Map) {
      return '{' + v.keySet().collect { it.toString() }.sort().collect { k -> "${k}:${canon(v[k], root)}" }.join(',') + '}'
    }
    if (v instanceof Collection) {
      def items = v.collect { canon(it, root) }
      // collect() delivers files in task-completion order; a set of files is the same input.
      if (!v.isEmpty() && v.every { it instanceof Path || it instanceof File }) items = items.sort()
      return '[' + items.join(',') + ']'
    }
    def s = v.toString()
    // A string that is an existing absolute file is a path that was passed as val.
    if (s.startsWith('/') && !s.contains('\n') && s.length() < 4096) {
      try {
        def p = Paths.get(s)
        if (Files.isRegularFile(p)) return fingerprint(p, root)
      } catch (Exception ignored) { }
    }
    return s.replace(root.toString(), '<OUT>')
  }

  static String fingerprint(Path p, Path root) {
    try {
      def real = p.toRealPath()
      def store = root.resolve('.store')
      def storeReal = Files.exists(store) ? store.toRealPath() : store
      if (real.startsWith(storeReal)) return 'store:' + storeReal.relativize(real)
      if (Files.isDirectory(real)) {
        def kids = []
        real.toFile().eachFileRecurse { f -> if (f.isFile()) kids << "${real.relativize(f.toPath())}:${f.length()}".toString() }
        return 'dir:' + real.fileName + ':' + sha(kids.sort().join(','))
      }
      long size = Files.size(real)
      if (size <= INLINE_MAX) return 'sha:' + sha(real.toFile().bytes)
      return "file:${real.fileName}:${size}".toString()
    } catch (Exception e) {
      return 'missing:' + p.fileName
    }
  }

  // ── helpers ──────────────────────────────────────────────────────────────

  static String sha(String s) { sha(s.getBytes('UTF-8')) }

  static String sha(byte[] b) {
    MessageDigest.getInstance('SHA-256').digest(b).encodeHex().toString()
  }

  private static void writeManifest(Path store, String module, String name, List<String> parts) {
    try {
      def f = store.resolve('_manifests').resolve(module).resolve("${name}.txt").toFile()
      if (!f.exists()) {
        f.parentFile.mkdirs()
        f.text = parts.join('\n') + '\n'
      }
    } catch (Exception ignored) { }   // manifest is a debugging aid only
  }
}
