// ============================================================================
// ResourceUtils.groovy — CPU/Memory helpers for TrackTx
// ============================================================================
//
// Moved here from nextflow.config because the strict syntax parser (default in
// Nextflow 26.04+) forbids function declarations inside config files.
// Reference these from nextflow.config as: ResourceUtils.alignmentCpus() etc.
// ============================================================================

class ResourceUtils {

    // ── System Detection ──────────────────────────────────────────────────────

    static def detectCpus() {
        def cpus = Runtime.runtime.availableProcessors()
        def envCpus = System.getenv('NXF_HOST_CPUS')
        return envCpus ? (envCpus as Integer) : cpus
    }

    static def detectMemoryGb() {
        try {
            def osBean = java.lang.management.ManagementFactory.operatingSystemMXBean
            def method = osBean.class.getMethod('getTotalPhysicalMemorySize')
            long totalBytes = (Long) method.invoke(osBean)
            return Math.max(2, (int)(totalBytes >> 30))
        } catch (Exception e) {
            def envMem = System.getenv('NXF_HOST_MEM')
            return envMem ? (envMem as Integer) : 8
        }
    }

    static def getSystemCpus()     { return Math.max(1, detectCpus()) }
    static def getSystemMemoryGb() { return Math.max(2, detectMemoryGb()) }

    // ── Tool-Specific CPU Allocation (based on efficiency benchmarks) ─────────

    static def alignmentCpus() {
        def total = getSystemCpus()
        if (total <= 4)  return total
        if (total <= 8)  return 4
        if (total <= 12) return 5
        if (total <= 16) return 6
        if (total <= 32) return 6
        return 8
    }

    static def alignmentForks() {
        def total = getSystemCpus()
        if (total < 8) return 1
        return Math.max(1, (int)(total / alignmentCpus()))
    }

    static def preprocessCpus() {
        def total = getSystemCpus()
        if (total <= 4) return Math.max(2, total - 1)
        return 4
    }

    static def preprocessForks() {
        def total = getSystemCpus()
        if (total < 6) return 1
        return Math.max(1, (int)(total / preprocessCpus()))
    }

    static def tracksCpus() {
        // 4 coverage jobs run in parallel (3'/5' ends × main/allMap BAMs); each is
        // single-threaded, so the task needs exactly 4 cores to keep all jobs busy.
        def total = getSystemCpus()
        if (total <= 4) return total
        return 4
    }

    static def tracksForks() {
        def total = getSystemCpus()
        if (total < 4) return 1
        return Math.max(1, (int)(total / tracksCpus()))
    }

    static def normCpus() {
        // 2 normalize_bedgraph jobs run in parallel per section; each is awk + bedGraphToBigWig
        // (both single-threaded), so 2 concurrent jobs is the ceiling.
        def total = getSystemCpus()
        if (total <= 4) return 2
        return 4
    }

    static def normForks() {
        def total = getSystemCpus()
        if (total < 10) return 1
        return Math.max(1, Math.min(3, (int)(total / normCpus())))
    }

    static def divergentCpus() {
        // Cores passed to --ncores; each sample parallelises across chromosomes.
        def total = getSystemCpus()
        if (total <= 4)  return total
        if (total <= 8)  return 4
        if (total <= 16) return 8
        if (total <= 32) return 8
        return 8
    }

    static def divergentForks() {
        def total = getSystemCpus()
        if (total < 8) return 1
        return Math.max(1, (int)(total / divergentCpus()))
    }

    static def indexCpus() {
        def total = getSystemCpus()
        if (total <= 4)  return total
        if (total <= 8)  return 6
        if (total <= 16) return 8
        if (total <= 32) return 12
        return 16
    }

    static def analysisCpus() {
        def total = getSystemCpus()
        if (total <= 4)  return 2
        if (total <= 8)  return 2
        if (total <= 16) return 3
        return 4
    }

    static def analysisForks() {
        def total = getSystemCpus()
        return Math.max(2, (int)(total / analysisCpus()))
    }

    static def lightweightCpus()  { return 1 }
    static def lightweightForks() { return Math.max(4, (int)(getSystemCpus() / 2)) }

    // ── Performance Profile Forks (external storage: USB, NAS, network) ──────

    static def performanceAlignmentForks() {
        def total = getSystemCpus()
        if (total < 8) return 2
        return Math.max(2, (int)(total / 4))
    }

    static def performancePreprocessForks() {
        def total = getSystemCpus()
        if (total < 6) return 2
        return Math.max(2, (int)(total / 3))
    }

    static def performanceTracksForks() {
        def total = getSystemCpus()
        if (total < 6) return 2
        return Math.max(2, (int)(total / 2))
    }

    static def performanceNormForks() {
        def total = getSystemCpus()
        if (total < 10) return 2
        return Math.max(2, Math.min(4, (int)(total / 4)))
    }

    static def performanceLightweightForks() {
        def total = getSystemCpus()
        return Math.max(6, (int)(total * 1.5))
    }

    static def performanceAnalysisForks() {
        def total = getSystemCpus()
        return Math.max(3, (int)(total / 2))
    }

    // ── Memory Allocation Strategies ─────────────────────────────────────────

    static def intensiveMemory(baseGb = 8) {
        def totalGb = getSystemMemoryGb()
        if (totalGb <= 8) {
            return "${Math.max(4, (int)(totalGb * 0.8))} GB"
        } else if (totalGb <= 32) {
            return "${Math.min((int)(baseGb * 2), (int)(totalGb * 0.5))} GB"
        } else {
            return "${Math.min((int)(baseGb * 3), (int)(totalGb * 0.6))} GB"
        }
    }

    static def moderateMemory(baseGb = 4) {
        def totalGb = getSystemMemoryGb()
        if (totalGb <= 8) {
            return "${Math.max(2, (int)(totalGb * 0.4))} GB"
        } else if (totalGb <= 32) {
            return "${Math.min((int)(baseGb * 1.5), (int)(totalGb * 0.3))} GB"
        } else {
            return "${Math.min((int)(baseGb * 2), (int)(totalGb * 0.4))} GB"
        }
    }

    static def lightMemory() {
        def totalGb = getSystemMemoryGb()
        if (totalGb <= 8)  return "1 GB"
        if (totalGb <= 32) return "2 GB"
        return "4 GB"
    }

    static def alignmentMemory() {
        def totalGb = getSystemMemoryGb()
        def forks = alignmentForks()
        if (totalGb <= 8) {
            return "${(int)(totalGb * 0.8)} GB"
        } else if (totalGb <= 32) {
            def perTask = Math.max(8, (int)((totalGb * 0.6) / Math.max(1, forks)))
            return "${Math.min(16, (int)perTask)} GB"
        } else {
            def perTask = Math.max(12, (int)((totalGb * 0.5) / Math.max(1, forks)))
            return "${Math.min(24, (int)perTask)} GB"
        }
    }
}
