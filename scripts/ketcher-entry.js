// ketcher-entry.js — source of the WebMMKetcher global inside ketcher-bundle.js.
// Reconstructed from the committed bundle (app/vendor/ketcher/ketcher-bundle.js):
// same API surface platform.html depends on — mount(el, staticResourcesUrl),
// unmount(), getMolfile(), setMolecule(mol). Do not rename these.
//
// The Indigo worker + wasm are embedded by ketcher-standalone's own dist
// (rollup web-worker-loader base64 blob); esbuild only needs to keep the
// string intact (it does — no special loader required).
import * as React from 'react';
import { createRoot } from 'react-dom/client';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';

let editorApi = null;   // Ketcher API once the editor signals onInit
let root = null;        // react-dom root, for unmount
let readyPromise = null;

window.WebMMKetcher = {
  mount(el, staticResourcesUrl) {
    const provider = new StandaloneStructServiceProvider();
    readyPromise = new Promise(resolve => {
      root = createRoot(el);
      root.render(React.createElement(Editor, {
        staticResourcesUrl: staticResourcesUrl || '',
        structServiceProvider: provider,
        onInit: api => { editorApi = api; resolve(api); }
      }));
    });
  },
  unmount() {
    if (root) root.unmount();
    root = null; editorApi = null; readyPromise = null;
  },
  async getMolfile() {
    return (await readyPromise).getMolfile();
  },
  async setMolecule(mol) {
    await (await readyPromise).setMolecule(mol);
  },
  // ---- Top-3 second half: query-as-drawing + selection bridge ----
  // The raw Ketcher instance, for anything not wrapped (feature-detect at
  // the call site — the exact surface varies across bundle versions).
  api() { return editorApi; },
  // SMARTS with the query features the user drew (any-atoms, ring bonds,
  // substitution counts). Throws when this bundle lacks it.
  async getSmarts() {
    const a = await readyPromise;
    if (!a || typeof a.getSmarts !== 'function') throw new Error('getSmarts unavailable in this ketcher build');
    return a.getSmarts();
  },
  // Current editor selection: { atoms?: number[], bonds?: number[] }
  // (1-based atom numbering in molfile order — feeds the atom ledger).
  async getSelection() {
    const a = await readyPromise;
    const ed = a && a.editor;
    if (!ed || typeof ed.selection !== 'function') throw new Error('selection unavailable in this ketcher build');
    const sel = ed.selection();
    if (!sel) return { atoms: [], bonds: [] };
    return { atoms: sel.atoms || [], bonds: sel.bonds || [] };
  }
};
