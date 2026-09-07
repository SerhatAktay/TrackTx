// Sanitize a user-supplied custom genome/spike-in id to safe path characters
// only, since it becomes part of a storeDir/cache path verbatim. Shared here
// (not duplicated per-file) so a future fix to the sanitization rule can't
// be applied in one call site and missed in another.
class GenomeId {
  static String sanitize(raw) {
    return raw.toString().trim().replaceAll(/[^A-Za-z0-9_.-]/, '_')
  }
}
