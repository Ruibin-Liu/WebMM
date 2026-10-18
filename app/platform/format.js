// platform/format.js — M4 P1 export/import format v1 (docs/project-management.md §2)
// FROZEN with format.version 1:
//   - canonJson: recursive key sort; numbers use JSON.stringify's own
//     toString semantics — NO toFixed rounding (review V3-1)
//   - volatile-field exclusion list is EMPTY as of v1 (commands are
//     deterministic; timestamps live in meta/facts, not commands)
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).Format = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function canonJson(v) {
    if (v === null || typeof v !== 'object') return JSON.stringify(v);
    if (Array.isArray(v)) return '[' + v.map(canonJson).join(',') + ']';
    return '{' + Object.keys(v).sort().map(k => JSON.stringify(k) + ':' + canonJson(v[k])).join(',') + '}';
  }

  function logHash(commands) {
    let h1 = 0x811c9dc5, h2 = 0x01000193;
    const str = canonJson(commands);
    for (let i = 0; i < str.length; i++) {
      h1 = Math.imul(h1 ^ str.charCodeAt(i), 16777619) >>> 0;
      h2 = (h2 + str.charCodeAt(i) * (i + 1)) >>> 0;
    }
    return ('0000000' + h1.toString(16)).slice(-8) + ('0000000' + h2.toString(16)).slice(-8);
  }

  // envelope builder — the page supplies raw {schemaVersion, meta, commands}
  // and layer-tagged facts; everything lineage-related is added here
  function buildExport(opts) {
    const raw = opts.raw, facts = opts.facts || [];
    const hash = logHash(raw.commands);
    return {
      format: 'webmm-project', version: 1,
      schemaVersion: raw.schemaVersion,
      sourceProjectId: opts.sourceProjectId || 'default',
      logHash: hash,
      projectMeta: { name: opts.name || '项目' },
      commands: raw.commands,
      facts,
      exportedAt: (opts.exportedAt || new Date()).toISOString ? (opts.exportedAt || new Date()).toISOString() : String(opts.exportedAt),
      engineVersions: opts.engineVersions || {},
    };
  }

  // parse + validate an import file; throws on structural problems.
  // Returns {commands, facts, schemaVersion, sourceProjectId, logHash,
  // hashOk} — hashOk is null when the file predates hashing (legacy).
  function parseImport(data) {
    if (!data || typeof data !== 'object') throw new Error('import: not a JSON object');
    const commands = Array.isArray(data.commands) ? data.commands : null;
    if (!commands) throw new Error('import: no commands array');
    if (data.format === 'webmm-project') {
      if (data.version !== 1) throw new Error('import: unsupported format version ' + data.version);
    } else if (data.format != null) {
      throw new Error('import: unknown format ' + data.format);
    }
    const facts = Array.isArray(data.facts) ? data.facts : [];
    const hashOk = data.logHash ? logHash(commands) === data.logHash : null;
    return {
      commands, facts,
      schemaVersion: data.schemaVersion || 1,
      sourceProjectId: data.sourceProjectId || null,
      logHash: data.logHash || null,
      hashOk,
      name: (data.projectMeta && data.projectMeta.name) || null,
    };
  }

  return { canonJson, logHash, buildExport, parseImport };
});
