# AI collaboration core — Capability / Plan / Tools / Evidence

> Status: first headless contract slice. No model vendor, UI, or browser wiring yet.
> Scope: one skill (`rgroup-optimization`) and one capability (`analogScan`).

## 1. Boundary

The model drafts and explains. WebMM validates, compiles, and executes. The model never reads DOM state, clicks UI, writes project storage, or calls WASM directly.

```text
model provider
  ↕ one planning call / one review call
Agent.Orchestrator
  ↕ allow-listed, sanitized tools
Capability registry + R-group skill + DesignPlan validator
  ↓ UI approval only
PlanCompiler → existing Tasks.Contract.makeSpec → existing task executor
  ↓ real spec + fact
EvidenceBundle → cited DesignReview → optional next plan draft
```

There is no autonomous loop. Execution happens only after the UI/execution layer issues an approval and supplies the current execution context.

## 2. Modules

| Module | Responsibility |
|---|---|
| `app/agent/capabilities.js` | Versioned, model-readable capability metadata, actual evidence types, budgets, and limitations |
| `app/agent/skills.js` | Versioned LBDD SOP for planning and reviewing R-group optimization |
| `app/agent/plans.js` | Plan normalization, scientific-goal checks, plan hash, approval, and TaskSpec compilation |
| `app/agent/tools.js` | Allow-listed model tools and model-boundary context sanitization |
| `app/agent/evidence.js` | Citation-ready artifacts from real task records and review-reference validation |
| `app/agent/orchestrator.js` | One bounded planning call and one bounded review call; no execution method |

All modules are UMD and dependency-free beyond the existing task contract. They run in browsers and Node.

## 3. Capability contract

`analogScan` revision 1 currently exposes only evidence produced by the real pipeline:

- descriptors already present on enumerated candidates;
- shape Tanimoto and optional flex Tanimoto;
- SA score when available.

It does **not** claim conformer-energy evidence, experimental activity prediction, binding affinity, solubility, or selectivity. Shape similarity is not activity evidence. MMFF energy is not a cross-molecule activity ranking.

`experimentalActivity`, `solubility`, and `selectivity` may appear in a plan only when also listed in `unverifiedGoals`. Unknown goals are rejected rather than silently relabelled.

## 4. DesignPlan v1

Required semantic fields:

```text
kind = rgroupOptimization
capability = analogScan
capabilityRevision = 1
skillId = rgroup-optimization
baseRevision / parentRef / siteRef / fragmentSet
hypothesis
goals / unverifiedGoals
budget { maxCandidates, refineTop, finalTop }
```

Unknown top-level and budget fields are rejected. Budgets are integers capped by the existing funnel at 400 / 30 / 10. Defaults shrink consistently for smaller runs.

`planHash` fingerprints the complete normalized plan using the existing task canonicalization. Editing goals, identity, hypothesis, or budget changes the hash.

## 5. Tools visible to the model

```text
getDesignContext()
listCapabilities()
getCapability(id)
getLibrarySummary(id)
validatePlan(plan)
getEvidence(taskId)
```

There is deliberately no `execute`, `approve`, or `compile` tool.

`getDesignContext` uses a whitelist. It may expose current design identity, parent SMILES, goals, audit, and candidate count. It removes parent SDF, candidate molblocks, build IDs, and other execution payloads. Review bundles contain result rows, not raw conformer SDF.

The host must treat remote-model sharing as a separate privacy decision. Sanitization prevents accidental over-sharing; it is not user consent.

## 6. Approval and compilation

Only the UI/execution layer calls:

```js
const approval = Agent.Plans.issueApproval(plan, { approvedBy: 'user:...' });
const spec = Agent.Plans.compileApprovedPlan({
  plan,
  approval,
  context: executionContext,
  taskId,
});
```

Approval binds:

- normalized `planHash`;
- `baseRevision`;
- non-empty `approvedBy`.

Compilation rechecks the plan against the current `revision`, `parentRef`, `siteRef`, and `fragmentSet`. It requires current parent structures, build ID, and enumerated candidates. Any edit, stale context, mismatched approval, missing identity, or excess candidate count aborts before task creation.

The compiler emits the existing frozen `analogScan input v1`. It changes only deterministic A/B/C/D caps. Conformer counts, seed offsets, ranking, and funnel behavior stay inherited from `Tasks.Contract`. Plan metadata, approval, and unverified goals live in `spec.provenance`, including `policyRev`.

## 7. Evidence and review

`Evidence.makeAnalogEvidence({ spec, fact })` copies real result values into:

- one summary `evidenceId`;
- one `evidenceId` per candidate;
- actual scores and descriptors;
- stage accounts and truncation;
- engine versions and policy revision;
- capability limitations and unverified goals.

`DesignReview.claims[]` requires non-empty text and at least one valid evidence reference. Unknown references are rejected. Counterexamples must reference candidate evidence, not the summary. `unresolved` preserves unmeasured or unevaluated questions.

A `nextPlanDraft` is validated against the current context and returned as `awaiting_approval`. Review never creates a user Decision or starts another run.

## 8. Provider interface

The orchestrator is model-vendor agnostic:

```ts
provider.plan({
  intent, context, capability, skill
}) -> DesignPlan JSON

provider.review({
  evidence, capability, skill
}) -> DesignReview JSON
```

A provider may call a local model, a remote model, or a deterministic test double. Network authentication, streaming, retries, and vendor tool protocols belong outside the core contract.

## 9. Tests and next slice

Run:

```sh
node tests/agent/run.js
node tests/platform/run.js
```

The next implementation slice should:

1. extract the existing R-group enumeration/context into a headless browser host;
2. load these modules in `platform.html` without changing current UI behavior;
3. connect approved compiled specs to the existing task store and funnel;
4. add a minimal Plan Card for edit/approve/run/review;
5. only then add a model-provider adapter.

No cloud provider, autonomous execution, general chat, additional skills, MCP, or new LBDD algorithms are part of this slice.
