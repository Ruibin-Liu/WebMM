// app/tasks/contract.js — Top-5 task contract (envelope layer, frozen in
// docs/task-contract.md). Pure: no I/O, no DOM, Node-testable.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Tasks = root.Tasks || {}).Contract = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const ENVELOPE_VERSION = 1;

  const KINDS = ['analogScan', 'conformers', 'optimize', 'shapeAlign'];
  const KIND_INPUT_FROZEN = { analogScan: 1 };   // others queue-only until frozen

  // state machine (§2) — every transition not listed is a contract error
  const TRANSITIONS = {
    queued: ['running', 'cancelled', 'abandoned'],
    running: ['complete', 'failed', 'cancelled', 'interrupted', 'abandoned'],
    complete: ['stale'],
    failed: [], cancelled: [], abandoned: [], stale: [], interrupted: [],
  };
  const STATES = Object.keys(TRANSITIONS);

  function canTransition(from, to) {
    return !!TRANSITIONS[from] && TRANSITIONS[from].includes(to);
  }

  // ---- canonical serialization (§4) ----
  const f64 = new Float64Array(1);
  const u64 = new BigUint64Array(f64.buffer);
  function numBits(n) {
    if (Number.isNaN(n)) return 'nan';
    if (n === Infinity) return '+inf';
    if (n === -Infinity) return '-inf';
    f64[0] = n;
    return u64[0].toString(16).padStart(16, '0');   // -0 and 0 differ in bits
  }
  // strings, numbers (IEEE754), booleans, null, arrays (order = semantics),
  // plain objects (sorted keys, undefined-valued keys dropped). Anything
  // else is a contract error, never silently coerced.
  function canon(v) {
    if (v === null) return 'null';
    const t = typeof v;
    if (t === 'number') return '#' + numBits(v);
    if (t === 'boolean') return v ? 'true' : 'false';
    if (t === 'string') return JSON.stringify(v);
    if (Array.isArray(v)) return '[' + v.map(canon).join(',') + ']';
    if (t === 'object') {
      const ks = Object.keys(v).filter(k => v[k] !== undefined).sort();
      return '{' + ks.map(k => JSON.stringify(k) + ':' + canon(v[k])).join(',') + '}';
    }
    throw new Error('canon: unsupported type ' + t + ' (functions/symbols are not fingerprintable)');
  }

  // FNV-1a 32 — same integer in every JS engine (UTF-16 code units)
  function fnv1a(str) {
    let h = 0x811c9dc5;
    for (let i = 0; i < str.length; i++) {
      h ^= str.charCodeAt(i);
      h = (h + ((h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24))) >>> 0;
    }
    return h >>> 0;
  }
  const hash = fnv1a;   // alias: contract speaks of fnv1a32 content hashes

  // ---- seed derivation with domain separation (§3) ----
  // Missing keys participate EXPLICITLY as '∅' — different call shapes must
  // never collide. stage must be a stable stageId/op name, never a position.
  function deriveSeed(parts) {
    const shape = {
      roundSeed: parts.roundSeed === undefined ? '∅' : parts.roundSeed,
      kind: parts.kind === undefined ? '∅' : parts.kind,
      stage: parts.stage === undefined ? '∅' : parts.stage,
      molId: parts.molId === undefined ? '∅' : parts.molId,
      confIdx: parts.confIdx === undefined ? '∅' : parts.confIdx,
    };
    return fnv1a(canon(shape));
  }

  // ---- effective input (§4): defaults expanded BEFORE fingerprinting ----
  const ANALOG_STAGES_V1 = [
    { id: 'A', op: 'enumerate', caps: { products: 400, rotbMax: 10, heavyDeltaMax: 9, heavyDeltaMin: -6 } },
    { id: 'B', op: 'score', nConfs: 10, seedOffset: 4200, cap: 400, inflightPerWorker: 2 },
    { id: 'C', op: 'score', nConfs: 30, seedOffset: 14200, cap: 30, combine: 'max' },
    { id: 'D', op: 'flex', cap: 10, skipSelfT: 0.9995 },
  ];

  function expandStages(stages) {
    if (!Array.isArray(stages) || !stages.length) return JSON.parse(JSON.stringify(ANALOG_STAGES_V1));
    // per-stage defaults over the declared stages (id/op required)
    return stages.map(s => {
      if (!s || typeof s.id !== 'string' || typeof s.op !== 'string') {
        throw new Error('expandStages: every stage needs string id and op');
      }
      const base = ANALOG_STAGES_V1.find(d => d.id === s.id && d.op === s.op) || {};
      const out = {};
      // declared fields first (deep copy of plain data), then fill defaults
      for (const k of Object.keys(s)) if (s[k] !== undefined) out[k] = JSON.parse(JSON.stringify(s[k]));
      for (const k of Object.keys(base)) if (out[k] === undefined) out[k] = JSON.parse(JSON.stringify(base[k]));
      return out;
    });
  }

  // minimal, project-specific validation (deliberately NOT a workflow framework)
  function validateStages(stages) {
    const OPS = { enumerate: 1, score: 1, flex: 1 };
    for (const s of stages) {
      if (!OPS[s.op]) throw new Error('validateStages: unknown op ' + s.op);
      if (s.cap !== undefined && (!Number.isInteger(s.cap) || s.cap < 1)) throw new Error('validateStages: cap must be a positive int');
      if (s.op === 'score' && (!Number.isInteger(s.nConfs) || s.nConfs < 1 || s.nConfs > 50)) throw new Error('validateStages: nConfs out of 1..50');
      if (s.op === 'score' && !Number.isInteger(s.seedOffset)) throw new Error('validateStages: seedOffset must be an int');
      if (s.op === 'score' && s.combine !== undefined && s.combine !== 'max') throw new Error('validateStages: combine must be "max"');
      if (s.op === 'flex' && !(s.skipSelfT >= 0 && s.skipSelfT <= 1)) throw new Error('validateStages: skipSelfT out of [0,1]');
    }
    // stage ids unique
    const ids = new Set(stages.map(s => s.id));
    if (ids.size !== stages.length) throw new Error('validateStages: duplicate stage ids');
    return true;
  }

  // effectiveInput per kind — the ONLY shape fingerprints and reruns see
  function effectiveInput(kind, input) {
    if (kind === 'analogScan') {
      const stages = expandStages(input && input.stages);
      validateStages(stages);
      return { ...input, stages };
    }
    return input;   // other kinds: queue-only, no expansion until frozen
  }

  function policyRev(stages) { return 'p' + fnv1a(canon(stages)).toString(16); }

  // fingerprint over {kind, inputVersion, effectiveInput, seed} (§4)
  function fingerprint(kind, inputVersion, effective, seed) {
    return 't' + fnv1a(canon({ kind, v: inputVersion, input: effective, seed: seed === undefined ? null : seed })).toString(16);
  }

  // ---- spec construction + guards ----
  function makeSpec(partial) {
    if (!KINDS.includes(partial.kind)) throw new Error('unknown kind ' + partial.kind);
    const eff = effectiveInput(partial.kind, partial.input || {});
    const inputVersion = KIND_INPUT_FROZEN[partial.kind] || 0;
    if (inputVersion === 0) throw new Error('kind ' + partial.kind + ' input schema not frozen yet (queue-only)');
    return {
      taskId: partial.taskId,
      kind: partial.kind,
      version: ENVELOPE_VERSION,
      input: eff,
      inputVersion,
      seed: partial.seed === undefined ? null : partial.seed,
      provenance: partial.provenance || {},
    };
  }

  function inputRevisionFor(kind, effInput) {
    if (kind === 'analogScan') {
      return {
        parent: 'h' + fnv1a(canon(effInput.parentSdf || '')).toString(16),
        fragmentSet: (effInput.inputRevision && effInput.inputRevision.fragmentSet) || 'builtin-v1',
        policy: policyRev(effInput.stages),
      };
    }
    return null;
  }

  return {
    ENVELOPE_VERSION, KINDS, KIND_INPUT_FROZEN, STATES, TRANSITIONS,
    canTransition, canon, fnv1a, hash, numBits, deriveSeed,
    ANALOG_STAGES_V1, expandStages, validateStages, effectiveInput,
    policyRev, fingerprint, makeSpec, inputRevisionFor,
  };
});
