let wasm;

function addToExternrefTable0(obj) {
    const idx = wasm.__externref_table_alloc();
    wasm.__wbindgen_externrefs.set(idx, obj);
    return idx;
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return decodeText(ptr, len);
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function handleError(f, args) {
    try {
        return f.apply(this, args);
    } catch (e) {
        const idx = addToExternrefTable0(e);
        wasm.__wbindgen_exn_store(idx);
    }
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    }
}

let WASM_VECTOR_LEN = 0;

const ConformerBatchFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_conformerbatch_free(ptr >>> 0, 1));

const ConvergenceOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_convergenceoptions_free(ptr >>> 0, 1));

const ETKDGResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_etkdgresult_free(ptr >>> 0, 1));

const MDLiveFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mdlive_free(ptr >>> 0, 1));

const MDOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mdoptions_free(ptr >>> 0, 1));

const MDResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mdresult_free(ptr >>> 0, 1));

const MetaDLiveFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_metadlive_free(ptr >>> 0, 1));

const MetaDOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_metadoptions_free(ptr >>> 0, 1));

const MetaDResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_metadresult_free(ptr >>> 0, 1));

const OptimizationOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_optimizationoptions_free(ptr >>> 0, 1));

const OptimizationResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_optimizationresult_free(ptr >>> 0, 1));

/**
 * Batch ETKDG embedding for the conformer ensemble: n structures with
 * seeds seed_base + i (diverse by construction; i=0 reproduces the
 * single-structure export's seed 42 when seed_base = 42).
 */
export class ConformerBatch {
    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(ConformerBatch.prototype);
        obj.__wbg_ptr = ptr;
        ConformerBatchFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ConformerBatchFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_conformerbatch_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get_n_atoms() {
        const ret = wasm.conformerbatch_get_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    get_n_confs() {
        const ret = wasm.conformerbatch_get_n_confs(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    get_success() {
        const ret = wasm.conformerbatch_get_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {Float64Array}
     */
    get_coordinates() {
        const ret = wasm.conformerbatch_get_coordinates(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {string}
     */
    get_error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.conformerbatch_get_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
}
if (Symbol.dispose) ConformerBatch.prototype[Symbol.dispose] = ConformerBatch.prototype.free;

/**
 * Convergence criteria options
 */
export class ConvergenceOptions {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ConvergenceOptionsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_convergenceoptions_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get max_force() {
        const ret = wasm.__wbg_get_convergenceoptions_max_force(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set max_force(arg0) {
        wasm.__wbg_set_convergenceoptions_max_force(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get rms_force() {
        const ret = wasm.__wbg_get_convergenceoptions_rms_force(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set rms_force(arg0) {
        wasm.__wbg_set_convergenceoptions_rms_force(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get energy_change() {
        const ret = wasm.__wbg_get_convergenceoptions_energy_change(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set energy_change(arg0) {
        wasm.__wbg_set_convergenceoptions_energy_change(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get max_iterations() {
        const ret = wasm.__wbg_get_convergenceoptions_max_iterations(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {number} arg0
     */
    set max_iterations(arg0) {
        wasm.__wbg_set_convergenceoptions_max_iterations(this.__wbg_ptr, arg0);
    }
}
if (Symbol.dispose) ConvergenceOptions.prototype[Symbol.dispose] = ConvergenceOptions.prototype.free;

export class ETKDGResult {
    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(ETKDGResult.prototype);
        obj.__wbg_ptr = ptr;
        ETKDGResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ETKDGResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_etkdgresult_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get_n_atoms() {
        const ret = wasm.conformerbatch_get_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    get_success() {
        const ret = wasm.etkdgresult_get_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {Float64Array}
     */
    get_coordinates() {
        const ret = wasm.etkdgresult_get_coordinates(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {string}
     */
    get_error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.etkdgresult_get_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
}
if (Symbol.dispose) ETKDGResult.prototype[Symbol.dispose] = ETKDGResult.prototype.free;

/**
 * Live MD handle: holds a running simulation so JS can advance it in small
 * chunks (one per animation frame) and render the trajectory as it evolves.
 */
export class MDLive {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MDLiveFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mdlive_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    steps_done() {
        const ret = wasm.mdlive_steps_done(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    temperature() {
        const ret = wasm.mdlive_temperature(this.__wbg_ptr);
        return ret;
    }
    /**
     * Per-atom force magnitudes (kcal/mol/Å) from the last force evaluation.
     * @returns {Float64Array}
     */
    force_magnitudes() {
        const ret = wasm.mdlive_force_magnitudes(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    potential_energy() {
        const ret = wasm.mdlive_potential_energy(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kinematically drag one atom: set its position, zero its velocity, and
     * refresh the cached energy/forces so the next step is consistent.
     * @param {number} i
     * @param {number} x
     * @param {number} y
     * @param {number} z
     */
    set_atom_position(i, x, y, z) {
        wasm.mdlive_set_atom_position(this.__wbg_ptr, i, x, y, z);
    }
    /**
     * Rescale velocities to a new target temperature (K) and retarget the
     * thermostat; re-initializes from Maxwell–Boltzmann if currently ~0 K;
     * `temperature_k <= 0` freezes the molecule.
     * @param {number} temperature_k
     */
    rescale_temperature(temperature_k) {
        wasm.mdlive_rescale_temperature(this.__wbg_ptr, temperature_k);
    }
    /**
     * @param {string} sdf_content
     * @param {MDOptions} options
     */
    constructor(sdf_content, options) {
        const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        _assertClass(options, MDOptions);
        var ptr1 = options.__destroy_into_raw();
        const ret = wasm.mdlive_new(ptr0, len0, ptr1);
        this.__wbg_ptr = ret >>> 0;
        MDLiveFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Advance the simulation by n steps (no-op if construction failed).
     * @param {number} n_steps
     */
    step(n_steps) {
        wasm.mdlive_step(this.__wbg_ptr, n_steps);
    }
    /**
     * @returns {string}
     */
    error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.mdlive_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Current coordinates as a flat [x0,y0,z0,x1,...] array.
     * @returns {Float64Array}
     */
    coords() {
        const ret = wasm.mdlive_coords(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Per-atom masses (amu) — for center-of-mass tracking in the demo.
     * @returns {Float64Array}
     */
    masses() {
        const ret = wasm.mdlive_masses(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    n_atoms() {
        const ret = wasm.mdlive_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    success() {
        const ret = wasm.mdlive_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    time_fs() {
        const ret = wasm.mdlive_time_fs(this.__wbg_ptr);
        return ret;
    }
}
if (Symbol.dispose) MDLive.prototype[Symbol.dispose] = MDLive.prototype.free;

/**
 * Options for running molecular dynamics (gas-phase MMFF).
 */
export class MDOptions {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MDOptionsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mdoptions_free(ptr, 0);
    }
    /**
     * @param {number} v
     */
    set_n_steps(v) {
        wasm.mdoptions_set_n_steps(this.__wbg_ptr, v);
    }
    /**
     * @param {string} v
     */
    set_mmff_variant(v) {
        const ptr0 = passStringToWasm0(v, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.mdoptions_set_mmff_variant(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {number} v
     */
    set_temperature_k(v) {
        wasm.mdoptions_set_temperature_k(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_friction_per_ps(v) {
        wasm.mdoptions_set_friction_per_ps(this.__wbg_ptr, v);
    }
    /**
     * Save a trajectory frame every N steps. 0 (default) = final frame only.
     * @param {number} v
     */
    set_snapshot_interval(v) {
        wasm.mdoptions_set_snapshot_interval(this.__wbg_ptr, v);
    }
    /**
     * Defaults: MMFF94s, 1 fs step, 1000 steps, 300 K, friction 1/ps, snapshots off, seed 42.
     */
    constructor() {
        const ret = wasm.mdoptions_new();
        this.__wbg_ptr = ret >>> 0;
        MDOptionsFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * @param {bigint} v
     */
    set_seed(v) {
        wasm.mdoptions_set_seed(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_dt_fs(v) {
        wasm.mdoptions_set_dt_fs(this.__wbg_ptr, v);
    }
}
if (Symbol.dispose) MDOptions.prototype[Symbol.dispose] = MDOptions.prototype.free;

/**
 * Result of an MD run: a sampled trajectory (flattened coordinates) plus per-frame
 * energies/temperatures and final stats.
 */
export class MDResult {
    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(MDResult.prototype);
        obj.__wbg_ptr = ptr;
        MDResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MDResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mdresult_free(ptr, 0);
    }
    /**
     * Flattened trajectory, n_frames * n_atoms * 3 (Float64Array in JS).
     * @returns {Float64Array}
     */
    coordinates() {
        const ret = wasm.mdresult_coordinates(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    final_energy() {
        const ret = wasm.mdresult_final_energy(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {Float64Array}
     */
    temperatures() {
        const ret = wasm.mdresult_temperatures(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    final_temperature() {
        const ret = wasm.mdresult_final_temperature(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {string}
     */
    error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.mdresult_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @returns {number}
     */
    steps() {
        const ret = wasm.mdresult_steps(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    n_atoms() {
        const ret = wasm.mdresult_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    success() {
        const ret = wasm.mdresult_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {Float64Array}
     */
    energies() {
        const ret = wasm.mdresult_energies(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    n_frames() {
        const ret = wasm.mdresult_n_frames(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {Float64Array}
     */
    times_fs() {
        const ret = wasm.mdresult_times_fs(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Single component: get_coord(frame, atom, axis), axis 0/1/2 = x/y/z.
     * @param {number} frame
     * @param {number} atom
     * @param {number} axis
     * @returns {number}
     */
    get_coord(frame, atom, axis) {
        const ret = wasm.mdresult_get_coord(this.__wbg_ptr, frame, atom, axis);
        return ret;
    }
}
if (Symbol.dispose) MDResult.prototype[Symbol.dispose] = MDResult.prototype.free;

/**
 * MMFF variant selection
 * @enum {0 | 1}
 */
export const MMFFVariant = Object.freeze({
    MMFF94: 0, "0": "MMFF94",
    MMFF94s: 1, "1": "MMFF94s",
});

/**
 * Live metadynamics handle: like [`MDLive`], plus CV/hill/FES accessors so
 * the demo can show the bias building up while the simulation runs.
 */
export class MetaDLive {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MetaDLiveFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_metadlive_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    hill_count() {
        const ret = wasm.metadlive_hill_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    steps_done() {
        const ret = wasm.mdlive_steps_done(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    temperature() {
        const ret = wasm.mdlive_temperature(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {Float64Array}
     */
    hill_centers() {
        const ret = wasm.metadlive_hill_centers(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    potential_energy() {
        const ret = wasm.mdlive_potential_energy(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {string} sdf_content
     * @param {MetaDOptions} options
     */
    constructor(sdf_content, options) {
        const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        _assertClass(options, MetaDOptions);
        var ptr1 = options.__destroy_into_raw();
        const ret = wasm.metadlive_new(ptr0, len0, ptr1);
        this.__wbg_ptr = ret >>> 0;
        MetaDLiveFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Advance the simulation by n steps (no-op if construction failed).
     * @param {number} n_steps
     */
    step(n_steps) {
        wasm.metadlive_step(this.__wbg_ptr, n_steps);
    }
    /**
     * @returns {string}
     */
    error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.metadlive_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @param {number} grid_points
     * @returns {Float64Array}
     */
    fes_f(grid_points) {
        const ret = wasm.metadlive_fes_f(this.__wbg_ptr, grid_points);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * FES grid points (CV values) and free energies for `grid_points` bins.
     * @param {number} grid_points
     * @returns {Float64Array}
     */
    fes_s(grid_points) {
        const ret = wasm.metadlive_fes_s(this.__wbg_ptr, grid_points);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    coords() {
        const ret = wasm.metadlive_coords(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Per-atom masses (amu) — for center-of-mass tracking in the demo.
     * @returns {Float64Array}
     */
    masses() {
        const ret = wasm.metadlive_masses(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * CV value of the last evaluated step.
     * @returns {number}
     */
    last_cv() {
        const ret = wasm.metadlive_last_cv(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {number}
     */
    n_atoms() {
        const ret = wasm.metadlive_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    success() {
        const ret = wasm.mdlive_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    time_fs() {
        const ret = wasm.mdlive_time_fs(this.__wbg_ptr);
        return ret;
    }
}
if (Symbol.dispose) MetaDLive.prototype[Symbol.dispose] = MetaDLive.prototype.free;

/**
 * Options for a metadynamics run (well-tempered, gas-phase MMFF).
 */
export class MetaDOptions {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MetaDOptionsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_metadoptions_free(ptr, 0);
    }
    /**
     * "dihedral" (4 atoms) or "distance" (2 atoms).
     * @param {string} v
     */
    set_cv_type(v) {
        const ptr0 = passStringToWasm0(v, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.metadoptions_set_cv_type(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {number} v
     */
    set_n_steps(v) {
        wasm.metadoptions_set_n_steps(this.__wbg_ptr, v);
    }
    /**
     * Comma-separated atom indices, e.g. "8,2,1,0" for a dihedral.
     * @param {string} v
     */
    set_cv_atoms(v) {
        const ptr0 = passStringToWasm0(v, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.metadoptions_set_cv_atoms(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {number} v
     */
    set_hill_width(v) {
        wasm.metadoptions_set_hill_width(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_bias_factor(v) {
        wasm.metadoptions_set_bias_factor(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_hill_height(v) {
        wasm.metadoptions_set_hill_height(this.__wbg_ptr, v);
    }
    /**
     * @param {string} v
     */
    set_mmff_variant(v) {
        const ptr0 = passStringToWasm0(v, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.metadoptions_set_mmff_variant(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {number} v
     */
    set_temperature_k(v) {
        wasm.mdoptions_set_temperature_k(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_fes_grid_points(v) {
        wasm.metadoptions_set_fes_grid_points(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_friction_per_ps(v) {
        wasm.mdoptions_set_friction_per_ps(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_deposit_interval(v) {
        wasm.metadoptions_set_deposit_interval(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_snapshot_interval(v) {
        wasm.metadoptions_set_snapshot_interval(this.__wbg_ptr, v);
    }
    constructor() {
        const ret = wasm.metadoptions_new();
        this.__wbg_ptr = ret >>> 0;
        MetaDOptionsFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * @param {number} v
     */
    set_seed(v) {
        wasm.metadoptions_set_seed(this.__wbg_ptr, v);
    }
    /**
     * @param {number} v
     */
    set_dt_fs(v) {
        wasm.mdoptions_set_dt_fs(this.__wbg_ptr, v);
    }
}
if (Symbol.dispose) MetaDOptions.prototype[Symbol.dispose] = MetaDOptions.prototype.free;

/**
 * Result of a metadynamics run: trajectory + CV trace + free-energy surface.
 */
export class MetaDResult {
    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(MetaDResult.prototype);
        obj.__wbg_ptr = ptr;
        MetaDResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MetaDResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_metadresult_free(ptr, 0);
    }
    /**
     * @returns {Float64Array}
     */
    coordinates() {
        const ret = wasm.metadresult_coordinates(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    final_energy() {
        const ret = wasm.mdresult_final_energy(this.__wbg_ptr);
        return ret;
    }
    /**
     * @returns {Float64Array}
     */
    hill_centers() {
        const ret = wasm.metadresult_hill_centers(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {string}
     */
    error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.metadresult_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Free energy (kcal/mol) at each grid point.
     * @returns {Float64Array}
     */
    fes_f() {
        const ret = wasm.metadresult_fes_f(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * CV grid points for the FES.
     * @returns {Float64Array}
     */
    fes_s() {
        const ret = wasm.metadresult_fes_s(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    n_atoms() {
        const ret = wasm.metadresult_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    n_hills() {
        const ret = wasm.metadresult_n_hills(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {boolean}
     */
    success() {
        const ret = wasm.metadresult_success(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {Float64Array}
     */
    energies() {
        const ret = wasm.metadresult_energies(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {number}
     */
    n_frames() {
        const ret = wasm.metadresult_n_frames(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {Float64Array}
     */
    times_fs() {
        const ret = wasm.metadresult_times_fs(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {Float64Array}
     */
    cv_values() {
        const ret = wasm.metadresult_cv_values(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
}
if (Symbol.dispose) MetaDResult.prototype[Symbol.dispose] = MetaDResult.prototype.free;

/**
 * Optimization options
 */
export class OptimizationOptions {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        OptimizationOptionsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_optimizationoptions_free(ptr, 0);
    }
    /**
     * @returns {string}
     */
    get mmff_variant() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.__wbg_get_optimizationoptions_mmff_variant(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @param {string} arg0
     */
    set mmff_variant(arg0) {
        const ptr0 = passStringToWasm0(arg0, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.__wbg_set_optimizationoptions_mmff_variant(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * Force-field engine: "MMFF94s" (default), "MMFF94", or "GFNFF".
     * Empty falls back to the legacy mmff_variant field.
     * @returns {string}
     */
    get engine() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.__wbg_get_optimizationoptions_engine(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Force-field engine: "MMFF94s" (default), "MMFF94", or "GFNFF".
     * Empty falls back to the legacy mmff_variant field.
     * @param {string} arg0
     */
    set engine(arg0) {
        const ptr0 = passStringToWasm0(arg0, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.__wbg_set_optimizationoptions_engine(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {number} val
     */
    set_max_force(val) {
        wasm.mdoptions_set_dt_fs(this.__wbg_ptr, val);
    }
    /**
     * @param {number} val
     */
    set_rms_force(val) {
        wasm.mdoptions_set_temperature_k(this.__wbg_ptr, val);
    }
    /**
     * @param {number} val
     */
    set_energy_change(val) {
        wasm.mdoptions_set_friction_per_ps(this.__wbg_ptr, val);
    }
    /**
     * @param {number} val
     */
    set_max_iterations(val) {
        wasm.optimizationoptions_set_max_iterations(this.__wbg_ptr, val);
    }
    constructor() {
        const ret = wasm.optimizationoptions_new();
        this.__wbg_ptr = ret >>> 0;
        OptimizationOptionsFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) OptimizationOptions.prototype[Symbol.dispose] = OptimizationOptions.prototype.free;

/**
 * Optimization result
 */
export class OptimizationResult {
    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(OptimizationResult.prototype);
        obj.__wbg_ptr = ptr;
        OptimizationResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        OptimizationResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_optimizationresult_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    get n_atoms() {
        const ret = wasm.__wbg_get_optimizationresult_n_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {number} arg0
     */
    set n_atoms(arg0) {
        wasm.__wbg_set_optimizationresult_n_atoms(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get final_energy() {
        const ret = wasm.__wbg_get_convergenceoptions_max_force(this.__wbg_ptr);
        return ret;
    }
    /**
     * @param {number} arg0
     */
    set final_energy(arg0) {
        wasm.__wbg_set_convergenceoptions_max_force(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {boolean}
     */
    get converged() {
        const ret = wasm.__wbg_get_optimizationresult_converged(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @param {boolean} arg0
     */
    set converged(arg0) {
        wasm.__wbg_set_optimizationresult_converged(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {number}
     */
    get iterations() {
        const ret = wasm.__wbg_get_optimizationresult_iterations(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {number} arg0
     */
    set iterations(arg0) {
        wasm.__wbg_set_optimizationresult_iterations(this.__wbg_ptr, arg0);
    }
    /**
     * @returns {string}
     */
    get message() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.__wbg_get_optimizationresult_message(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @param {string} arg0
     */
    set message(arg0) {
        const ptr0 = passStringToWasm0(arg0, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.__wbg_set_optimizationresult_message(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * Engine actually used: "MMFF94s", "MMFF94", or "GFNFF".
     * @returns {string}
     */
    get_engine() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.optimizationresult_get_engine(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Per-atom partial charges at the optimized geometry.
     * @returns {Float64Array}
     */
    get_charges() {
        const ret = wasm.optimizationresult_get_charges(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @returns {string}
     */
    get_message() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.optimizationresult_get_message(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @returns {boolean}
     */
    get_converged() {
        const ret = wasm.optimizationresult_get_converged(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * @returns {number}
     */
    get_iterations() {
        const ret = wasm.optimizationresult_get_iterations(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {number}
     */
    get_final_energy() {
        const ret = wasm.mdresult_final_energy(this.__wbg_ptr);
        return ret;
    }
    /**
     * JSON object of per-term energies (kcal/mol), engine-specific keys.
     * @returns {string}
     */
    get_energy_terms_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.optimizationresult_get_energy_terms_json(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @param {number} atom_idx
     * @param {number} coord_idx
     * @returns {number}
     */
    get_coord(atom_idx, coord_idx) {
        const ret = wasm.optimizationresult_get_coord(this.__wbg_ptr, atom_idx, coord_idx);
        return ret;
    }
}
if (Symbol.dispose) OptimizationResult.prototype[Symbol.dispose] = OptimizationResult.prototype.free;

/**
 * Add explicit hydrogens to a heavy-atom SDF/molblock (v1.1): H counts from
 * standard organic valences (aromatic = 1.5, charge-adjusted), single bonds
 * to the parent heavy atom, coordinates placed away from the neighbor
 * centroid (L-BFGS relaxes them; ETKDG regenerates everything on embed).
 * Inputs that already contain explicit H are returned unchanged.
 * @param {string} sdf_content
 * @returns {string}
 */
export function add_hydrogens_wasm(sdf_content) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.add_hydrogens_wasm(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Attach hydrogens to a 3D heavy-atom SDF with geometry-aware placement
 * (tetrahedral fans / bisectors at standard bond lengths). Pairs with the
 * heavy-only ETKDG embed, which is ~20x faster than embedding all-H.
 * @param {string} sdf_content
 * @returns {string}
 */
export function attach_hydrogens_3d_wasm(sdf_content) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.attach_hydrogens_3d_wasm(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Single-point energy + term decomposition for a 3D structure (no
 * optimization). Returns a JSON string: {"E", "terms": {...}, "engine"}.
 * Used by the workbench to sync the energy panel when switching conformers.
 * @param {string} sdf_content
 * @param {string} engine
 * @returns {string}
 */
export function energy_terms_wasm(sdf_content, engine) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(engine, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.energy_terms_wasm(ptr0, len0, ptr1, len1);
        var ptr3 = ret[0];
        var len3 = ret[1];
        if (ret[3]) {
            ptr3 = 0; len3 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
    }
}

/**
 * @param {string} sdf_content
 * @param {number} n
 * @param {bigint} seed_base
 * @returns {ConformerBatch}
 */
export function generate_conformers_wasm(sdf_content, n, seed_base) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.generate_conformers_wasm(ptr0, len0, n, seed_base);
    return ConformerBatch.__wrap(ret);
}

/**
 * @param {string} sdf_content
 * @returns {ETKDGResult}
 */
export function generate_initial_coordinates_wasm(sdf_content) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.generate_initial_coordinates_wasm(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return ETKDGResult.__wrap(ret[0]);
}

export function init() {
    wasm.init();
}

/**
 * @param {string} sdf_content
 * @param {OptimizationOptions} options
 * @returns {OptimizationResult}
 */
export function optimize_from_sdf(sdf_content, options) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    _assertClass(options, OptimizationOptions);
    var ptr1 = options.__destroy_into_raw();
    const ret = wasm.optimize_from_sdf(ptr0, len0, ptr1);
    return OptimizationResult.__wrap(ret);
}

/**
 * Optimize from an SDF using the SDF coordinates directly (no ETKDG).
 * If the SDF is 2D, coordinates are used as-is (z=0 plane).
 * @param {string} sdf_content
 * @param {OptimizationOptions} options
 * @returns {OptimizationResult}
 */
export function optimize_from_sdf_direct(sdf_content, options) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    _assertClass(options, OptimizationOptions);
    var ptr1 = options.__destroy_into_raw();
    const ret = wasm.optimize_from_sdf_direct(ptr0, len0, ptr1);
    return OptimizationResult.__wrap(ret);
}

/**
 * Run molecular dynamics on an SDF molecule (gas-phase MMFF) and return a sampled
 * trajectory. NVT (BAOAB Langevin) if `friction_per_ps > 0`, else NVE.
 * @param {string} sdf_content
 * @param {MDOptions} options
 * @returns {MDResult}
 */
export function run_md_from_sdf(sdf_content, options) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    _assertClass(options, MDOptions);
    var ptr1 = options.__destroy_into_raw();
    const ret = wasm.run_md_from_sdf(ptr0, len0, ptr1);
    return MDResult.__wrap(ret);
}

/**
 * Run well-tempered metadynamics on an SDF molecule (gas-phase MMFF) and return
 * the trajectory + CV trace + free-energy surface (FES).
 * @param {string} sdf_content
 * @param {MetaDOptions} options
 * @returns {MetaDResult}
 */
export function run_metadynamics_from_sdf(sdf_content, options) {
    const ptr0 = passStringToWasm0(sdf_content, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    _assertClass(options, MetaDOptions);
    var ptr1 = options.__destroy_into_raw();
    const ret = wasm.run_metadynamics_from_sdf(ptr0, len0, ptr1);
    return MetaDResult.__wrap(ret);
}

/**
 * WebMM crate version (CARGO_PKG_VERSION), exposed for export provenance.
 * @returns {string}
 */
export function webmm_version() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.webmm_version();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

const EXPECTED_RESPONSE_TYPES = new Set(['basic', 'cors', 'default']);

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && EXPECTED_RESPONSE_TYPES.has(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else {
                    throw e;
                }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }
}

function __wbg_get_imports() {
    const imports = {};
    imports.wbg = {};
    imports.wbg.__wbg___wbindgen_is_undefined_f6b95eab589e0269 = function(arg0) {
        const ret = arg0 === undefined;
        return ret;
    };
    imports.wbg.__wbg___wbindgen_throw_dd24417ed36fc46e = function(arg0, arg1) {
        throw new Error(getStringFromWasm0(arg0, arg1));
    };
    imports.wbg.__wbg_call_abb4ff46ce38be40 = function() { return handleError(function (arg0, arg1) {
        const ret = arg0.call(arg1);
        return ret;
    }, arguments) };
    imports.wbg.__wbg_error_7534b8e9a36f1ab4 = function(arg0, arg1) {
        let deferred0_0;
        let deferred0_1;
        try {
            deferred0_0 = arg0;
            deferred0_1 = arg1;
            console.error(getStringFromWasm0(arg0, arg1));
        } finally {
            wasm.__wbindgen_free(deferred0_0, deferred0_1, 1);
        }
    };
    imports.wbg.__wbg_new_8a6f238a6ece86ea = function() {
        const ret = new Error();
        return ret;
    };
    imports.wbg.__wbg_new_no_args_cb138f77cf6151ee = function(arg0, arg1) {
        const ret = new Function(getStringFromWasm0(arg0, arg1));
        return ret;
    };
    imports.wbg.__wbg_now_2c95c9de01293173 = function(arg0) {
        const ret = arg0.now();
        return ret;
    };
    imports.wbg.__wbg_performance_7a3ffd0b17f663ad = function(arg0) {
        const ret = arg0.performance;
        return ret;
    };
    imports.wbg.__wbg_stack_0ed75d68575b0f3c = function(arg0, arg1) {
        const ret = arg1.stack;
        const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
        getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
    };
    imports.wbg.__wbg_static_accessor_GLOBAL_769e6b65d6557335 = function() {
        const ret = typeof global === 'undefined' ? null : global;
        return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
    };
    imports.wbg.__wbg_static_accessor_GLOBAL_THIS_60cf02db4de8e1c1 = function() {
        const ret = typeof globalThis === 'undefined' ? null : globalThis;
        return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
    };
    imports.wbg.__wbg_static_accessor_SELF_08f5a74c69739274 = function() {
        const ret = typeof self === 'undefined' ? null : self;
        return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
    };
    imports.wbg.__wbg_static_accessor_WINDOW_a8924b26aa92d024 = function() {
        const ret = typeof window === 'undefined' ? null : window;
        return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
    };
    imports.wbg.__wbindgen_cast_2241b6af4c4b2941 = function(arg0, arg1) {
        // Cast intrinsic for `Ref(String) -> Externref`.
        const ret = getStringFromWasm0(arg0, arg1);
        return ret;
    };
    imports.wbg.__wbindgen_init_externref_table = function() {
        const table = wasm.__wbindgen_externrefs;
        const offset = table.grow(4);
        table.set(0, undefined);
        table.set(offset + 0, undefined);
        table.set(offset + 1, null);
        table.set(offset + 2, true);
        table.set(offset + 3, false);
    };

    return imports;
}

function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    __wbg_init.__wbindgen_wasm_module = module;
    cachedDataViewMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;


    wasm.__wbindgen_start();
    return wasm;
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (typeof module !== 'undefined') {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (typeof module_or_path !== 'undefined') {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (typeof module_or_path === 'undefined') {
        module_or_path = new URL('webmm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync };
export default __wbg_init;
