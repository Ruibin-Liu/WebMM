#!/usr/bin/env node
// Local-only WebMM Agent development server: static files + pi/magpie proxy.
// No request/response bodies are logged (molecular structures may be private).
import http from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(fileURLToPath(new URL('..', import.meta.url)));
const PORT = Number(process.env.AGENT_DEV_PORT || 8902);
const MAGPIE_URL = process.env.MAGPIE_URL || 'http://127.0.0.1:3425';
const MODEL = process.env.AGENT_MODEL || 'zcode/GLM-5.3-Flash';
const MAX_BODY = 1024 * 1024;
const MIME = {
  '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8', '.json': 'application/json; charset=utf-8', '.svg': 'image/svg+xml',
  '.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg', '.gif': 'image/gif',
  '.wasm': 'application/wasm', '.bin': 'application/octet-stream', '.txt': 'text/plain; charset=utf-8',
};

function json(res, status, value) {
  const body = JSON.stringify(value);
  res.writeHead(status, { 'Content-Type': 'application/json; charset=utf-8', 'Content-Length': Buffer.byteLength(body), 'Cache-Control': 'no-store' });
  res.end(body);
}

async function readJsonBody(req) {
  let size = 0;
  const chunks = [];
  for await (const chunk of req) {
    size += chunk.length;
    if (size > MAX_BODY) { const error = new Error('request body exceeds 1 MiB'); error.code = 'BODY_TOO_LARGE'; throw error; }
    chunks.push(chunk);
  }
  try { return JSON.parse(Buffer.concat(chunks).toString('utf8')); }
  catch (error) { error.code = 'INVALID_JSON'; throw error; }
}

async function agentApi(req, res, pathname) {
  if (pathname === '/api/agent/status' && req.method === 'GET') {
    json(res, 200, { ok: true, provider: 'pi/magpie', providerRevision: 'magpie-provider-v1', model: MODEL, contextRoute: 'pi/magpie' });
    return true;
  }
  if (pathname === '/api/agent/chat' && req.method === 'POST') {
    try {
      const body = await readJsonBody(req);
      if (!Array.isArray(body.messages) || !body.messages.length || body.messages.length > 12) {
        json(res, 400, { error: { code: 'MESSAGES_INVALID', message: 'messages must contain 1..12 entries' } });
        return true;
      }
      for (const message of body.messages) {
        if (!message || (message.role !== 'system' && message.role !== 'user' && message.role !== 'assistant') || typeof message.content !== 'string') {
          json(res, 400, { error: { code: 'MESSAGES_INVALID', message: 'each message needs a valid role and string content' } });
          return true;
        }
      }
      const upstreamBody = {
        model: MODEL,
        messages: body.messages,
        temperature: Number.isFinite(body.temperature) ? body.temperature : 0.1,
        max_tokens: Number.isInteger(body.max_tokens) ? Math.min(Math.max(body.max_tokens, 1), 8192) : 4096,
        stream: false,
      };
      const headers = { 'Content-Type': 'application/json' };
      if (process.env.MAGPIE_API_KEY) headers.Authorization = 'Bearer ' + process.env.MAGPIE_API_KEY;
      const upstream = await fetch(MAGPIE_URL + '/v1/chat/completions', {
        method: 'POST', headers, body: JSON.stringify(upstreamBody), signal: AbortSignal.timeout(180000),
      });
      const text = await upstream.text();
      res.writeHead(upstream.status, { 'Content-Type': 'application/json; charset=utf-8', 'Cache-Control': 'no-store' });
      res.end(text);
    } catch (error) {
      const status = error.code === 'BODY_TOO_LARGE' ? 413 : error.code === 'INVALID_JSON' ? 400 : 502;
      json(res, status, { error: { code: error.code || 'UPSTREAM_ERROR', message: error.message } });
    }
    return true;
  }
  if (pathname.startsWith('/api/agent/')) { json(res, 404, { error: { code: 'NOT_FOUND', message: 'unknown agent endpoint' } }); return true; }
  return false;
}

async function serveStatic(res, pathname) {
  let decoded;
  try { decoded = decodeURIComponent(pathname); }
  catch { json(res, 400, { error: { code: 'BAD_PATH', message: 'invalid URL encoding' } }); return; }
  const relative = decoded === '/' ? 'app/platform.html' : decoded.replace(/^\/+/, '');
  if (relative.split('/').some(part => part === '..' || part === '.git' || part === 'node_modules' || part === 'target')) {
    json(res, 403, { error: { code: 'FORBIDDEN', message: 'path is not served' } });
    return;
  }
  const filename = path.resolve(ROOT, relative);
  if (!filename.startsWith(ROOT + path.sep)) { json(res, 403, { error: { code: 'FORBIDDEN', message: 'path is not served' } }); return; }
  try {
    const info = await stat(filename);
    if (!info.isFile()) throw new Error('not a file');
    const data = await readFile(filename);
    res.writeHead(200, {
      'Content-Type': MIME[path.extname(filename)] || 'application/octet-stream',
      'Content-Length': data.length,
      'Cache-Control': 'no-store',
      'X-Content-Type-Options': 'nosniff',
    });
    res.end(data);
  } catch { json(res, 404, { error: { code: 'NOT_FOUND', message: 'file not found' } }); }
}

const server = http.createServer(async (req, res) => {
  const pathname = new URL(req.url, 'http://localhost').pathname;
  if (await agentApi(req, res, pathname)) return;
  if (req.method !== 'GET' && req.method !== 'HEAD') { json(res, 405, { error: { code: 'METHOD', message: 'method not allowed' } }); return; }
  await serveStatic(res, pathname);
});
server.listen(PORT, '127.0.0.1', () => {
  console.log(`WebMM Agent: http://127.0.0.1:${PORT}/app/platform.html`);
  console.log(`Model route: pi/magpie ${MODEL} (payloads are not logged)`);
});
for (const signal of ['SIGINT', 'SIGTERM']) process.on(signal, () => server.close(() => process.exit(0)));
