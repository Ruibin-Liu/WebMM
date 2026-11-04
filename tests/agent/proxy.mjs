// Local agent proxy tests: real HTTP proxy, fake magpie upstream.
import http from 'node:http';
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(fileURLToPath(new URL('../..', import.meta.url)));
let pass = 0, fail = 0;
function check(name, condition, detail) {
  if (condition) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + String(detail).slice(0, 160) : '')); }
}
const listen = (server, port) => new Promise((resolve, reject) => { server.once('error', reject); server.listen(port, '127.0.0.1', resolve); });
const close = server => new Promise(resolve => server.close(resolve));

let upstream, child;
try {
  const requests = [];
  upstream = http.createServer(async (req, res) => {
    const chunks = [];
    for await (const chunk of req) chunks.push(chunk);
    requests.push({ url: req.url, body: JSON.parse(Buffer.concat(chunks).toString('utf8')) });
    res.writeHead(200, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify({ choices: [{ message: { content: '{"ok":true}' } }] }));
  });
  await listen(upstream, 0);
  const upstreamPort = upstream.address().port;
  const devPort = 18000 + process.pid % 1000;
  child = spawn(process.execPath, ['scripts/agent-dev-server.mjs'], {
    cwd: ROOT,
    env: { ...process.env, AGENT_DEV_PORT: String(devPort), MAGPIE_URL: `http://127.0.0.1:${upstreamPort}`, AGENT_MODEL: 'zcode/Test-GLM' },
    stdio: ['ignore', 'pipe', 'pipe'],
  });
  let stderr = '';
  child.stderr.on('data', chunk => { stderr += chunk; });
  const base = `http://127.0.0.1:${devPort}`;
  let status;
  for (let i = 0; i < 60; i++) {
    try { status = await fetch(base + '/api/agent/status').then(r => r.json()); break; }
    catch { await new Promise(resolve => setTimeout(resolve, 100)); }
  }
  check('status exposes the server-selected model', status && status.model === 'zcode/Test-GLM', stderr);

  const chat = await fetch(base + '/api/agent/chat', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ model: 'attacker/model', messages: [{ role: 'user', content: 'test' }], max_tokens: 999999 }),
  });
  const chatData = await chat.json();
  check('chat forwards only a successful OpenAI-shaped response', chat.status === 200 && chatData.choices[0].message.content === '{"ok":true}');
  check('browser cannot override model or unbounded token budget',
    requests.length === 1 && requests[0].body.model === 'zcode/Test-GLM' && requests[0].body.max_tokens === 8192 && requests[0].body.stream === false,
    JSON.stringify(requests[0] && requests[0].body));
  check('message content reaches upstream unchanged', requests[0].body.messages[0].content === 'test');

  const traversal = await fetch(base + '/..%2f..%2fetc%2fpasswd');
  check('encoded path traversal is rejected', traversal.status === 403, traversal.status);
  const unknown = await fetch(base + '/api/agent/unknown');
  check('unknown agent endpoints return JSON 404', unknown.status === 404 && (await unknown.json()).error.code === 'NOT_FOUND');
  const invalid = await fetch(base + '/api/agent/chat', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{bad' });
  check('invalid JSON is rejected before upstream', invalid.status === 400 && requests.length === 1);
  const tooMany = await fetch(base + '/api/agent/chat', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ messages: Array.from({ length: 13 }, () => ({ role: 'user', content: 'x' })) }),
  });
  check('message count is bounded', tooMany.status === 400 && requests.length === 1);
  const page = await fetch(base + '/app/platform.html');
  check('static server serves the Agent-enabled platform', page.status === 200 && (await page.text()).includes('magpie-provider.js'));
} catch (error) {
  fail++; console.error(error);
} finally {
  if (child) child.kill('SIGTERM');
  if (upstream) await close(upstream);
}
console.log(`RESULT: ${pass} passed, ${fail} failed`);
process.exitCode = fail ? 1 : 0;
