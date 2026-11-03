# AI collaboration core — Capability / Plan / Tools / Evidence

> Status: browser vertical slice complete with a mock provider. No real model vendor is configured.
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
| `app/agent/rgroup.js` | Headless stage-A enumeration shared by the legacy explorer and Agent |
| `app/agent/browser-host.js` | Current R-group identity, execution context, enumeration, library summary, and real evidence access |
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
baseRevision / parentRef / siteRef / fragmentSet / fragmentCategory
hypothesis
goals / unverifiedGoals
budget { maxCandidates, refineTop, finalTop }
```

Unknown top-level and budget fields are rejected. `fragmentCategory` defaults to `all` and is part of both `planHash` and the browser context revision; changing the category after approval invalidates that approval. Budgets are integers capped by the existing funnel at 400 / 30 / 10. Defaults shrink consistently for smaller runs.

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

## 8. Browser execution and Plan Card

`platform.html` loads the Agent modules and creates `WebMMAgentHost` / `WebMMAgentTools`. The legacy explorer and approved Agent run both call `Agent.RGroup.enumerateRgroup`, then share the same worker pool, scheduler, task store, funnel, and Done-fact transaction.

The minimal Plan Card supports:

1. local draft, or an injected `window.WebMMAgentProvider`;
2. hypothesis, budget, and unverified-activity editing;
3. local approval followed by headless enumeration and `runAnalogFromCands`;
4. local evidence summary, or provider review with real evidence references.

Without a provider, the UI explicitly says “local draft” and “local evidence summary (not AI review)”. It never simulates an AI response. With a provider, planning and review each make one bounded call. Provider text reaches the page only through form values or `textContent`.

## 9. Provider interface

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

## 10. GLM feasibility probe

`glm-5.3-flash` was tested through Sider with the R-group contract:

- its first DesignPlan was schema-valid, including an explicit unverified experimental-activity goal;
- its first DesignReview added unrequested fields and unsupported interpretation, so the deterministic validator rejected it;
- after a schema-specific corrective prompt, it produced a valid review with only real evidence references and no unsupported activity claim.

GLM is therefore a candidate provider for **structured draft + deterministic validation + corrective retry**, not a trusted execution authority. The Sider session is only a test harness. Production use still needs a separately configured GLM-compatible endpoint or a thin model proxy and an explicit privacy decision.

## 11. Tests and next slice

Run:

```sh
node tests/agent/run.js
node tests/platform/run.js
tests/cdp/run.sh agent_plan
```

`agent_plan` uses a schema-valid mock provider and verifies the complete browser path plus the legacy explorer. `m12_tasks` remains the determinism, cancellation, and reconciliation regression gate.

The next implementation slice should:

1. add a configurable GLM/OpenAI-compatible provider with schema-specific corrective retry;
2. define endpoint/key handling and exactly which molecular context may leave the browser;
3. persist Plan/Review as first-class project artifacts rather than only task provenance;
4. add a second LBDD skill only after the R-group loop is used successfully.

No autonomous execution, general chat, MCP, or new LBDD algorithms are part of this slice.
