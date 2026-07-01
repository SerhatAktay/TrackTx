// ============================================================================
// DEPRECATED — intentionally empty.
// ============================================================================
//
// The CPU/memory helpers that used to live here have been moved into the
// `params.res` closure map in nextflow.config.
//
// Reason: under the Nextflow 26 strict config parser (ConfigDsl v2), classes
// from this lib/ directory are NOT visible to config closures. Referencing
// `ResourceUtils.x()` from a process resource directive therefore failed at
// task-submission time with:
//
//     No such variable: ResourceUtils
//
// Do not re-add resource helpers here. Edit `params.res` in nextflow.config.
// ============================================================================
