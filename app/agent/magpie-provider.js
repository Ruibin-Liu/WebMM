// app/agent/magpie-provider.js — real pi/magpie provider through the local
// relative proxy. Schema validation and one corrective retry happen here;
// execution remains exclusively with the approval-bound local layer.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./plans.js'), require('./evidence.js'));
  else { (root.Agent = root.Agent || {}).MagpieProvider = factory(root.Agent.Plans, root.Agent.Evidence); }
})(typeof self !== 'undefined' ? self : this, function (Plans, Evidence) {
  'use strict';

  function planRules(request) {
    const capability = request.capability, skill = request.skill, context = request.context;
    return `Return ONLY one JSON object with exactly these keys: version, kind, capability, capabilityRevision, skillId, skillVersion, baseRevision, parentRef, siteRef, fragmentSet, fragmentCategory, hypothesis, goals, unverifiedGoals, budget. Use version=1, kind=${capability.planKind}, capability=${capability.id}, capabilityRevision=${capability.revision}, skillId=${skill.id}, skillVersion=${skill.version}. Copy baseRevision, parentRef, siteRef, fragmentSet=${context.fragmentSet}, fragmentCategory=${context.fragmentCategory} from the supplied context. budget has only maxCandidates, refineTop, finalTop with limits 400/30/10 and maxCandidates>=refineTop>=finalTop. Computable goals are only ${capability.supportedGoals.join(', ')}. ${capability.knownUnverifiedGoals.join(', ')} must also appear in unverifiedGoals. Never claim predicted activity. Respect these capability limitations: ${capability.limitations.join(' ')}`;
  }
  const REVIEW_RULES = `Return ONLY one JSON object with exactly these keys: claims, counterexamples, unresolved, nextPlanDraft. claims is [{text:string,evidenceRefs:string[]}], and every claim needs at least one real evidenceId from the supplied evidence. counterexamples is a string array containing only candidate evidenceIds. unresolved is a string array. nextPlanDraft must be null. Do not add fields, invent values or references, infer molecular identity, or claim experimental activity. Treat truncated=true as incomplete ranking.`;

  function textOf(message) {
    const content = message && message.content;
    if (typeof content === 'string') return content;
    if (Array.isArray(content)) return content.map(part => typeof part === 'string' ? part : part && (part.text || part.content) || '').join('');
    return '';
  }

  function parseJsonText(text) {
    const trimmed = text.trim();
    const unfenced = trimmed.replace(/^```(?:json)?\s*/i, '').replace(/\s*```$/, '').trim();
    try { return { value: JSON.parse(unfenced), raw: unfenced }; }
    catch { /* recover a JSON object embedded in prose */ }
    const start = unfenced.indexOf('{');
    if (start < 0) throw new Error('model response contains no JSON object');
    let depth = 0, quoted = false, escaped = false;
    for (let i = start; i < unfenced.length; i++) {
      const ch = unfenced[i];
      if (escaped) { escaped = false; continue; }
      if (ch === '\\' && quoted) { escaped = true; continue; }
      if (ch === '"') { quoted = !quoted; continue; }
      if (quoted) continue;
      if (ch === '{') depth++;
      else if (ch === '}' && --depth === 0) {
        const raw = unfenced.slice(start, i + 1);
        return { value: JSON.parse(raw), raw };
      }
    }
    throw new Error('model response has an unterminated JSON object');
  }

  function createMagpieProvider(options) {
    const endpoint = options && options.endpoint || '/api/agent';
    const fetchImpl = options && options.fetchImpl || (typeof fetch !== 'undefined' ? fetch.bind(globalThis) : null);
    if (!fetchImpl) throw new Error('fetch is required');

    async function request(pathname, init) {
      const response = await fetchImpl(endpoint + pathname, init);
      const text = await response.text();
      let data;
      try { data = JSON.parse(text); }
      catch { throw new Error('agent endpoint returned non-JSON status ' + response.status + ': ' + text.slice(0, 160)); }
      if (!response.ok) throw new Error(data && data.error && data.error.message || 'agent endpoint status ' + response.status);
      return data;
    }

    async function status() { return request('/status', { method: 'GET' }); }

    async function chat(messages, maxTokens) {
      const data = await request('/chat', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ messages, temperature: 0.1, max_tokens: maxTokens }),
      });
      const choice = data && data.choices && data.choices[0];
      const text = textOf(choice && choice.message);
      if (!text.trim()) {
        const usage = data && data.usage && data.usage.completion_tokens;
        throw new Error('model returned empty content (finish_reason=' + (choice && choice.finish_reason || 'unknown') + ', completion_tokens=' + (usage == null ? 'unknown' : usage) + ')');
      }
      return parseJsonText(text);
    }

    async function structured(messages, validate, maxTokens) {
      let current = messages.slice();
      let lastErrors = [];
      for (let attempt = 0; attempt < 2; attempt++) {
        const parsed = await chat(current, maxTokens);
        const validation = validate(parsed.value);
        if (validation.ok) return validation;
        lastErrors = validation.errors;
        if (attempt === 0) {
          current = current.concat([
            { role: 'assistant', content: parsed.raw },
            { role: 'user', content: 'The JSON failed deterministic validation with these errors: ' + lastErrors.join('; ') + '. Return only a corrected JSON object that satisfies every rule and keeps all evidence references real.' },
          ]);
        }
      }
      const error = new Error('model JSON failed validation after one corrective retry: ' + lastErrors.join('; '));
      error.code = 'PROVIDER_SCHEMA';
      throw error;
    }

    async function plan(request) {
      const validation = await structured([
        { role: 'system', content: planRules(request) },
        { role: 'user', content: JSON.stringify({ intent: request.intent, context: request.context, capability: request.capability, skill: request.skill }) },
      ], raw => Plans.validatePlan(raw, request.context), 4096);
      return validation.plan;
    }

    async function review(request) {
      const validation = await structured([
        { role: 'system', content: REVIEW_RULES },
        { role: 'user', content: JSON.stringify({ evidence: request.evidence, capability: request.capability, skill: request.skill }) },
      ], raw => Evidence.validateReview(raw, request.evidence), 8192);
      return validation.review;
    }

    return Object.freeze({ status, plan, review });
  }

  return { createMagpieProvider, parseJsonText };
});
