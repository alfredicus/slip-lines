/**
 * Stress Field Computation Module
 * 
 * Implements theoretical slip lines and equipotentials from the Python implementation.
 * Based on: Taboada, A., 1993, Stress and strain from striated pebbles.
 */

// ============================================================================
// TYPES
// ============================================================================

export interface Point3D {
    x: number;
    y: number;
    z: number;
}

export interface SlipLineCurve {
    x: number[];
    y: number[];
    z: number[];
    normal_stress: number[];
    shear_stress: number[];
    friction: number[];
    type: 'slip_line' | 'meridian';
}

export interface EquipotentialCurve {
    x: number[];
    y: number[];
    z: number[];
    iso_level: number;
    normal_stress: number[];
    type: 'equipotential' | 'sigma2_meridian';
}

export interface PrincipalAxis {
    axis: Point3D;
    value: number;
    colorHex: number;
}

export interface ComputeConfig {
    nCurves?: number;
    nIsoCurves?: number;
    curveResolution?: number;
}

// ============================================================================
// CONSTANTS
// ============================================================================

const REVOLUTION_TOLERANCE = 1e-6;
const DEFAULT_N_CURVES = 12;
const DEFAULT_N_ISO = 8;
const DEFAULT_CURVE_RESOLUTION = 200;

export function getEquipotentialPoint(
    at: Point3D,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
): number {
    const x = at.x;
    const y = at.y;
    const z = at.z;

    // Convert Cartesian to spherical
    const r = Math.sqrt(x * x + y * y + z * z);
    const theta_rad = Math.acos(Math.max(-1, Math.min(1, z / r)));
    let phi_rad = Math.atan2(y, x);
    if (phi_rad < 0) phi_rad += 2 * Math.PI;

    const theta_deg = (theta_rad * 180) / Math.PI;
    const phi_deg = (phi_rad * 180) / Math.PI;

    // Calculate normal stress at this point
    return calculateNormalStress(phi_deg, theta_deg, lambda_x, lambda_y, lambda_z);
}

// ============================================================================
// COORDINATE TRANSFORMATIONS
// ============================================================================

/**
 * Convert spherical coordinates (phi, theta) to Cartesian (x, y, z)
 * 
 * phi: azimuth angle in degrees [0, 360]
 * theta: colatitude angle in degrees [0, 180]
 * r: radius (default 1.0)
 */
export function sphericalToCartesian(
    phi_deg: number,
    theta_deg: number,
    r: number = 1.0
): Point3D {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;

    return {
        x: r * Math.sin(theta_rad) * Math.cos(phi_rad),
        y: r * Math.sin(theta_rad) * Math.sin(phi_rad),
        z: r * Math.cos(theta_rad)
    };
}

/**
 * Convert multiple spherical coordinates to Cartesian (vectorized)
 */
function sphericalToCartesianVec(
    phi_deg: number[],
    theta_deg: number[]
): { x: number[]; y: number[]; z: number[] } {
    const x: number[] = [];
    const y: number[] = [];
    const z: number[] = [];

    for (let i = 0; i < phi_deg.length; i++) {
        const pt = sphericalToCartesian(phi_deg[i], theta_deg[i]);
        x.push(pt.x);
        y.push(pt.y);
        z.push(pt.z);
    }

    return { x, y, z };
}

// ============================================================================
// STRESS CALCULATIONS
// ============================================================================

/**
 * Calculate normal stress magnitude on sphere surface
 * 
 * Formula: Fn = (λx*cos²φ + λy*sin²φ)*sin²θ + λz*cos²θ
 */
export function calculateNormalStress(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;

    const cos_phi = Math.cos(phi_rad);
    const sin_phi = Math.sin(phi_rad);
    const cos_theta = Math.cos(theta_rad);
    const sin_theta = Math.sin(theta_rad);

    return (
        (lambda_x * cos_phi * cos_phi + lambda_y * sin_phi * sin_phi) * sin_theta * sin_theta +
        lambda_z * cos_theta * cos_theta
    );
}

/**
 * Calculate normal stress magnitude (vectorized)
 */
function calculateNormalStressVec(
    phi_deg: number[],
    theta_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number[] {
    const result: number[] = [];
    for (let i = 0; i < phi_deg.length; i++) {
        result.push(calculateNormalStress(phi_deg[i], theta_deg[i], lambda_x, lambda_y, lambda_z));
    }
    return result;
}

/**
 * Calculate shear stress magnitude on sphere surface
 */
function calculateShearStress(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;

    const cos_phi = Math.cos(phi_rad);
    const sin_phi = Math.sin(phi_rad);
    const cos_theta = Math.cos(theta_rad);
    const sin_theta = Math.sin(theta_rad);

    // Unit normal vector
    const nx = sin_theta * cos_phi;
    const ny = sin_theta * sin_phi;
    const nz = cos_theta;

    // Stress vector: σ·n
    const sx = lambda_x * nx;
    const sy = lambda_y * ny;
    const sz = lambda_z * nz;

    // Normal stress component
    const normal = lambda_x * nx * nx + lambda_y * ny * ny + lambda_z * nz * nz;

    // Shear vector = stress vector - normal component
    const shear_x = sx - normal * nx;
    const shear_y = sy - normal * ny;
    const shear_z = sz - normal * nz;

    return Math.sqrt(shear_x * shear_x + shear_y * shear_y + shear_z * shear_z);
}

/**
 * Calculate shear stress magnitude (vectorized)
 */
function calculateShearStressVec(
    phi_deg: number[],
    theta_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number[] {
    const result: number[] = [];
    for (let i = 0; i < phi_deg.length; i++) {
        result.push(calculateShearStress(phi_deg[i], theta_deg[i], lambda_x, lambda_y, lambda_z));
    }
    return result;
}

/**
 * Calculate friction coefficient (shear/normal stress ratio)
 */
function calculateFriction(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const normal = calculateNormalStress(phi_deg, theta_deg, lambda_x, lambda_y, lambda_z);
    const shear = calculateShearStress(phi_deg, theta_deg, lambda_x, lambda_y, lambda_z);

    if (Math.abs(normal) < 1e-10) return 0;
    return shear / Math.abs(normal);
}

/**
 * Calculate friction coefficient (vectorized)
 */
function calculateFrictionVec(
    phi_deg: number[],
    theta_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number[] {
    const result: number[] = [];
    for (let i = 0; i < phi_deg.length; i++) {
        result.push(calculateFriction(phi_deg[i], theta_deg[i], lambda_x, lambda_y, lambda_z));
    }
    return result;
}

// ============================================================================
// SLIP LINE EQUATION
// ============================================================================

/**
 * Slip line equation: tan(θ) = k1 * (sin φ)^exp_sin * (cos φ)^exp_cos
 * 
 * Returns theta in degrees given phi in degrees
 */
export function slipLineEquation(
    phi_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k1: number = 1.0
): number {
    // Check for revolution case
    if (Math.abs(lambda_y - lambda_x) < REVOLUTION_TOLERANCE) {
        return 45.0; // Special case: meridians
    }

    const phi_rad = (phi_deg * Math.PI) / 180;

    // Stress ratios
    const exp_sin = (lambda_x - lambda_z) / (lambda_y - lambda_x);
    const exp_cos = (lambda_z - lambda_y) / (lambda_y - lambda_x);

    // Trigonometric components with numerical stability
    const sin_phi = Math.max(Math.abs(Math.sin(phi_rad)), 1e-10);
    const cos_phi = Math.max(Math.abs(Math.cos(phi_rad)), 1e-10);

    // Slip line equation
    const tan_theta = k1 * Math.pow(sin_phi, exp_sin) * Math.pow(cos_phi, exp_cos);
    const theta_rad = Math.atan(tan_theta);
    const theta_deg = (theta_rad * 180) / Math.PI;

    return theta_deg;
}

/**
 * Slip line equation (vectorized)
 */
function slipLineEquationVec(
    phi_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k1: number = 1.0
): number[] {
    const result: number[] = [];
    for (let i = 0; i < phi_deg.length; i++) {
        result.push(slipLineEquation(phi_deg[i], lambda_x, lambda_y, lambda_z, k1));
    }
    return result;
}

// ============================================================================
// CURVE REFLECTION
// ============================================================================

/**
 * Reflect curve from one octant to all 8 octants using symmetry
 */
function reflectCurveToAllOctants(
    phi_deg: number[],
    theta_deg: number[]
): Array<{ phi: number[]; theta: number[]; x: number[]; y: number[]; z: number[] }> {
    const { x: x_base, y: y_base, z: z_base } = sphericalToCartesianVec(phi_deg, theta_deg);

    const reflections = [
        { sx: 1, sy: 1, sz: 1 },   // Original
        { sx: -1, sy: 1, sz: 1 },  // Reflect X
        { sx: 1, sy: -1, sz: 1 },  // Reflect Y
        { sx: 1, sy: 1, sz: -1 },  // Reflect Z
        { sx: -1, sy: -1, sz: 1 }, // Reflect X, Y
        { sx: -1, sy: 1, sz: -1 }, // Reflect X, Z
        { sx: 1, sy: -1, sz: -1 }, // Reflect Y, Z
        { sx: -1, sy: -1, sz: -1 } // Reflect X, Y, Z
    ];

    const curves: Array<{ phi: number[]; theta: number[]; x: number[]; y: number[]; z: number[] }> = [];

    for (const { sx, sy, sz } of reflections) {
        const x_refl = x_base.map((v) => sx * v);
        const y_refl = y_base.map((v) => sy * v);
        const z_refl = z_base.map((v) => sz * v);

        const phi_refl: number[] = [];
        const theta_refl: number[] = [];

        for (let i = 0; i < x_refl.length; i++) {
            const r = Math.sqrt(x_refl[i] * x_refl[i] + y_refl[i] * y_refl[i] + z_refl[i] * z_refl[i]);
            theta_refl.push((Math.acos(Math.max(-1, Math.min(1, z_refl[i] / r))) * 180) / Math.PI);

            let phi = (Math.atan2(y_refl[i], x_refl[i]) * 180) / Math.PI;
            if (phi < 0) phi += 360;
            phi_refl.push(phi);
        }

        curves.push({
            phi: phi_refl,
            theta: theta_refl,
            x: x_refl,
            y: y_refl,
            z: z_refl
        });
    }

    return curves;
}

// ============================================================================
// SLIP LINE GENERATION
// ============================================================================

/**
 * Generate slip lines for the sphere
 * Main implementation from Python's generate_slip_lines_3d
 */
export function generateSlipLines(
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    config?: Partial<ComputeConfig>
): SlipLineCurve[] {
    const cfg = {
        nCurves: config?.nCurves ?? DEFAULT_N_CURVES,
        curveResolution: config?.curveResolution ?? DEFAULT_CURVE_RESOLUTION
    };

    const allCurves: SlipLineCurve[] = [];

    // Create phi ranges with higher density at edges
    const createPhiRange = (n: number): number[] => {
        const phi: number[] = [];
        // Dense near poles: [0.5, 10]
        for (let i = 0; i < 50; i++) {
            phi.push(0.5 + (i / 50) * 9.5);
        }
        // Medium in middle: [10, 80]
        for (let i = 0; i < n - 100; i++) {
            phi.push(10 + (i / (n - 100)) * 70);
        }
        // Dense near pole: [80, 89.5]
        for (let i = 0; i < 50; i++) {
            phi.push(80 + (i / 50) * 9.5);
        }
        return phi;
    };

    // Check for revolution case
    if (Math.abs(lambda_y - lambda_x) < REVOLUTION_TOLERANCE) {
        // Revolution case: generate meridians
        for (let i = 0; i < cfg.nCurves; i++) {
            const azimuth = 0.5 + (i * (90 - 0.5)) / cfg.nCurves;
            const theta_range = createPhiRange(cfg.curveResolution);
            const phi_range = theta_range.map(() => azimuth);

            const { x, y, z } = sphericalToCartesianVec(phi_range, theta_range);
            const normal_stress = calculateNormalStressVec(phi_range, theta_range, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_range, theta_range, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_range, theta_range, lambda_x, lambda_y, lambda_z);

            // Reflect to all octants
            const reflected = reflectCurveToAllOctants(phi_range, theta_range);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'meridian'
                });
            }
        }
    } else {
        // General case: use k1 values with slip line equation
        const phi_range = createPhiRange(cfg.curveResolution);

        // Generate k1 values logarithmically distributed
        const k1_values: number[] = [];
        for (let i = 0; i < cfg.nCurves; i++) {
            const logVal = Math.log10(0.1) + (i / (cfg.nCurves - 1)) * (Math.log10(10) - Math.log10(0.1));
            k1_values.push(Math.pow(10, logVal));
        }

        for (const k1 of k1_values) {
            const theta_range = slipLineEquationVec(phi_range, lambda_x, lambda_y, lambda_z, k1);

            // Filter valid points
            const valid: boolean[] = [];
            for (let i = 0; i < theta_range.length; i++) {
                valid.push(theta_range[i] >= 0.1 && theta_range[i] <= 89.9);
            }

            const validCount = valid.filter((v) => v).length;
            if (validCount < 5) continue;

            // Extract valid points
            const phi_valid: number[] = [];
            const theta_valid: number[] = [];
            for (let i = 0; i < valid.length; i++) {
                if (valid[i]) {
                    phi_valid.push(phi_range[i]);
                    theta_valid.push(theta_range[i]);
                }
            }

            const { x, y, z } = sphericalToCartesianVec(phi_valid, theta_valid);
            const normal_stress = calculateNormalStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);

            // Reflect to all octants
            const reflected = reflectCurveToAllOctants(phi_valid, theta_valid);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'slip_line'
                });
            }
        }
    }

    return allCurves;
}

// ============================================================================
// EQUIPOTENTIAL GENERATION
// ============================================================================


/**
 * Generate equipotential curves (iso-normal stress)
 */
export function generateEquipotentials(
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    config?: Partial<ComputeConfig>
): EquipotentialCurve[] {
    const cfg = {
        nIsoCurves: config?.nIsoCurves ?? DEFAULT_N_ISO,
        curveResolution: config?.curveResolution ?? DEFAULT_CURVE_RESOLUTION
    };

    const allCurves: EquipotentialCurve[] = [];

    // Sort stress values
    const stresses = [lambda_x, lambda_y, lambda_z].sort((a, b) => a - b);
    const sigma1 = stresses[0]; // minimum (most compressive)
    const sigma3 = stresses[2]; // maximum (most extensional)

    // Generate iso levels with linear spacing
    const fn_levels: number[] = [];
    for (let i = 0; i <= cfg.nIsoCurves; i++) {
        const fn = sigma3 + ((sigma1 - sigma3) * i) / cfg.nIsoCurves;
        fn_levels.push(fn);
    }

    // Generate equipotential curves
    for (const target_fn of fn_levels) {
        // Skip endpoints (singular points)
        if (Math.abs(target_fn - sigma3) < 1e-6 || Math.abs(target_fn - sigma1) < 1e-6) {
            continue;
        }

        // Generate by sweeping through theta and finding phi values
        const theta_range: number[] = [];
        for (let i = 0; i < cfg.curveResolution; i++) {
            theta_range.push((i / (cfg.curveResolution - 1)) * 179.8 + 0.1);
        }

        const phi_values: number[] = [];
        const valid_theta: number[] = [];
        const valid_phi: number[] = [];

        for (const theta of theta_range) {
            // Binary search for phi where stress equals target_fn
            let phi_low = 0;
            let phi_high = 90;

            for (let iter = 0; iter < 30; iter++) {
                const phi_mid = (phi_low + phi_high) / 2;
                const stress = calculateNormalStress(phi_mid, theta, lambda_x, lambda_y, lambda_z);

                if (stress < target_fn) {
                    phi_low = phi_mid;
                } else {
                    phi_high = phi_mid;
                }
            }

            const phi = (phi_low + phi_high) / 2;
            const final_stress = calculateNormalStress(phi, theta, lambda_x, lambda_y, lambda_z);

            // Accept if close enough to target
            if (Math.abs(final_stress - target_fn) < Math.abs(sigma1 - sigma3) * 1e-3) {
                valid_theta.push(theta);
                valid_phi.push(phi);
            }
        }

        if (valid_theta.length > 5) {
            const { x, y, z } = sphericalToCartesianVec(valid_phi, valid_theta);
            const normal_stress = calculateNormalStressVec(valid_phi, valid_theta, lambda_x, lambda_y, lambda_z);

            allCurves.push({
                x,
                y,
                z,
                iso_level: target_fn,
                normal_stress,
                type: 'equipotential'
            });
        }
    }

    return allCurves;
}

// ============================================================================
// PRINCIPAL AXES
// ============================================================================

/**
 * Generate principal stress axes
 */
export function generatePrincipalAxes(
    sigma_x: number,
    sigma_y: number,
    sigma_z: number
): PrincipalAxis[] {
    const stresses = [
        { value: sigma_x, axis: { x: 1, y: 0, z: 0 } },
        { value: sigma_y, axis: { x: 0, y: 1, z: 0 } },
        { value: sigma_z, axis: { x: 0, y: 0, z: 1 } }
    ];

    // Sort by magnitude (descending)
    stresses.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));

    const colors = [0xff0000, 0x00ff00, 0x0000ff]; // Red, Green, Blue

    return stresses.map((stress, i) => ({
        axis: stress.axis,
        value: stress.value,
        colorHex: colors[i]
    }));
}