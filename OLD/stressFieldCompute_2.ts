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
    isSpecial?: boolean;
    specialType?: 'great_circle_xy' | 'great_circle_xz' | 'great_circle_yz';
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
 * Formula: Fn = (Î»x*cosÂ²Ï† + Î»y*sinÂ²Ï†)*sinÂ²Î¸ + Î»z*cosÂ²Î¸
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

    // Stress vector: ÏƒÂ·n
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
 * Slip line equation: tan(Î¸) = k1 * (sin Ï†)^exp_sin * (cos Ï†)^exp_cos
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
// HELPER FUNCTIONS FOR SLIP LINE GENERATION
// ============================================================================

/**
 * Solve for the integration constant k at a specific point on the slip line (general case)
 * tan(θ) = k × sin(φ)^exp_sin × cos(φ)^exp_cos
 */
function solveForKFromPoint(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;
    
    const exp_sin = (lambda_x - lambda_z) / (lambda_y - lambda_x);
    const exp_cos = (lambda_z - lambda_y) / (lambda_y - lambda_x);
    
    const sin_phi = Math.max(Math.abs(Math.sin(phi_rad)), 1e-10);
    const cos_phi = Math.max(Math.abs(Math.cos(phi_rad)), 1e-10);
    const tan_theta = Math.tan(theta_rad);
    
    const denominator = Math.pow(sin_phi, exp_sin) * Math.pow(cos_phi, exp_cos);
    const k = Math.abs(tan_theta) / Math.max(denominator, 1e-10);
    
    return k;
}

/**
 * Solve for k at a point on R=0 meridian using simplified equation
 * R=0: exp_cos=0, so cos_phi^0=1 (avoid undefined)
 * Simplified: tan(θ) = k × sin(φ)^exp_sin
 */
function solveForKFromPointR0(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;
    
    const exp_sin = (lambda_x - lambda_z) / (lambda_y - lambda_x); // = -1 for R=0
    
    const sin_phi = Math.max(Math.abs(Math.sin(phi_rad)), 1e-10);
    const tan_theta = Math.tan(theta_rad);
    
    // Only sin_phi term: k = tan(θ) / sin(φ)^exp_sin
    const denominator = Math.pow(sin_phi, exp_sin);
    const k = Math.abs(tan_theta) / Math.max(denominator, 1e-10);
    
    return k;
}

/**
 * Solve for k at a point on R=1 meridian using simplified equation
 * R=1: exp_sin=0, so sin_phi^0=1 (avoid undefined)
 * Simplified: tan(θ) = k × cos(φ)^exp_cos
 */
function solveForKFromPointR1(
    phi_deg: number,
    theta_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;
    
    const exp_cos = (lambda_z - lambda_y) / (lambda_y - lambda_x); // = -1 for R=1
    
    const cos_phi = Math.max(Math.abs(Math.cos(phi_rad)), 1e-10);
    const tan_theta = Math.tan(theta_rad);
    
    // Only cos_phi term: k = tan(θ) / cos(φ)^exp_cos
    const denominator = Math.pow(cos_phi, exp_cos);
    const k = Math.abs(tan_theta) / Math.max(denominator, 1e-10);
    
    return k;
}

/**
 * Compute slip line equation for R=0 (only sin_phi term)
 * Simplified: tan(θ) = k × sin(φ)^exp_sin
 */
function slipLineEquationR0(
    phi_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k: number = 1.0
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const exp_sin = (lambda_x - lambda_z) / (lambda_y - lambda_x);
    
    const sin_phi = Math.max(Math.abs(Math.sin(phi_rad)), 1e-10);
    const tan_theta = k * Math.pow(sin_phi, exp_sin);
    const theta_rad = Math.atan(tan_theta);
    
    return (theta_rad * 180) / Math.PI;
}

/**
 * Compute slip line equation for R=1 (only cos_phi term)
 * Simplified: tan(θ) = k × cos(φ)^exp_cos
 */
function slipLineEquationR1(
    phi_deg: number,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k: number = 1.0
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const exp_cos = (lambda_z - lambda_y) / (lambda_y - lambda_x);
    
    const cos_phi = Math.max(Math.abs(Math.cos(phi_rad)), 1e-10);
    const tan_theta = k * Math.pow(cos_phi, exp_cos);
    const theta_rad = Math.atan(tan_theta);
    
    return (theta_rad * 180) / Math.PI;
}

/**
 * Vectorized version of slipLineEquationR0
 */
function slipLineEquationVecR0(
    phi_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k: number = 1.0
): number[] {
    return phi_deg.map(phi => slipLineEquationR0(phi, lambda_x, lambda_y, lambda_z, k));
}

/**
 * Vectorized version of slipLineEquationR1
 */
function slipLineEquationVecR1(
    phi_deg: number[],
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    k: number = 1.0
): number[] {
    return phi_deg.map(phi => slipLineEquationR1(phi, lambda_x, lambda_y, lambda_z, k));
}

// ============================================================================
// SLIP LINE GENERATION
// ============================================================================


/**
 * Generate slip lines for the sphere
 * Handles three cases: R=0 (σ2=σ3), R=1 (σ2=σ1), and general case
 * Uses equidistant meridian points with solved k values for homogeneous distribution
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

    // CASE 1: R = 0 (σ2 = σ3, i.e., lambda_z ≈ lambda_y)
    // Revolution around σ1 (x-axis), meridians at φ=90° (y,z plane)
    if (Math.abs(lambda_z - lambda_y) < REVOLUTION_TOLERANCE) {
        console.log('[Slip Lines] R=0 case: σ2 ≈ σ3');
        console.log(`  lambda_z=${lambda_z}, lambda_y=${lambda_y}, diff=${Math.abs(lambda_z - lambda_y)}`);
        
        // Generate equidistant theta points on (y,z) meridian
        const theta_meridian: number[] = [];
        for (let i = 1; i < cfg.nCurves; i++) {
            const t = i / cfg.nCurves;
            theta_meridian.push(t * 90); // θ from 0° to 90°
        }
        console.log(`  Generated ${theta_meridian.length} equidistant meridian points:`, theta_meridian.map(t => t.toFixed(1)));
        
        // Fixed azimuth for revolution case: φ = 90° (y,z plane)
        const phi_fixed = 90;
        
        // Solve for k at each equidistant theta point using R=0 simplified equation
        for (const theta_point of theta_meridian) {
            const k = solveForKFromPointR0(phi_fixed, theta_point, lambda_x, lambda_y, lambda_z);
            console.log(`  At θ=${theta_point.toFixed(1)}°, φ=${phi_fixed}°: k=${k.toFixed(6)}`);
            
            // Generate slip line with this k value using R=0 simplified equation
            const phi_vals = [];
            for (let j = 0; j < cfg.curveResolution; j++) {
                phi_vals.push(0.1 + (j / (cfg.curveResolution - 1)) * 89.8);
            }
            
            // Use simplified R=0 equation (only sin_phi term)
            const theta_vals = slipLineEquationVecR0(phi_vals, lambda_x, lambda_y, lambda_z, k);
            const valid: boolean[] = theta_vals.map(t => t >= 0.1 && t <= 89.9);
            const validCount = valid.filter(v => v).length;
            console.log(`    Generated ${theta_vals.length} theta values, ${validCount} valid (θ in [0.1, 89.9])°`);
            if (validCount < 5) {
                console.log(`    Skipping: too few valid points`);
                continue;
            }

            const phi_valid: number[] = [];
            const theta_valid: number[] = [];
            for (let i = 0; i < valid.length; i++) {
                if (valid[i]) {
                    phi_valid.push(phi_vals[i]);
                    theta_valid.push(theta_vals[i]);
                }
            }

            const { x, y, z } = sphericalToCartesianVec(phi_valid, theta_valid);
            const normal_stress = calculateNormalStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);

            const reflected = reflectCurveToAllOctants(phi_valid, theta_valid);
            console.log(`    Created ${reflected.length} reflected curves (8 octants)`);
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

        // Add special slip lines: great circles (x,y) and (x,z)
        const specialCurves = generateSpecialSlipLinesR0(lambda_x, lambda_y, lambda_z, cfg.curveResolution);
        console.log(`  Added ${specialCurves.length} special slip lines (great circles)`);
        allCurves.push(...specialCurves);

    // CASE 2: R = 1 (σ2 = σ1, i.e., lambda_z ≈ lambda_x)
    // Revolution around σ3 (y-axis), meridians at φ=0° (x,z plane)
    } else if (Math.abs(lambda_z - lambda_x) < REVOLUTION_TOLERANCE) {
        console.log('[Slip Lines] R=1 case: σ2 ≈ σ1');
        console.log(`  lambda_z=${lambda_z}, lambda_x=${lambda_x}, diff=${Math.abs(lambda_z - lambda_x)}`);
        
        // Generate equidistant theta points on (x,z) meridian
        const theta_meridian: number[] = [];
        for (let i = 1; i < cfg.nCurves; i++) {
            const t = i / cfg.nCurves;
            theta_meridian.push(t * 90); // θ from 0° to 90°
        }
        console.log(`  Generated ${theta_meridian.length} equidistant meridian points:`, theta_meridian.map(t => t.toFixed(1)));
        
        // Fixed azimuth for revolution case: φ = 0° (x,z plane)
        const phi_fixed = 0;
        
        // Solve for k at each equidistant theta point using R=1 simplified equation
        for (const theta_point of theta_meridian) {
            const k = solveForKFromPointR1(phi_fixed, theta_point, lambda_x, lambda_y, lambda_z);
            console.log(`  At θ=${theta_point.toFixed(1)}°, φ=${phi_fixed}°: k=${k.toFixed(6)}`);
            
            // Generate slip line with this k value using R=1 simplified equation
            const phi_vals = [];
            for (let j = 0; j < cfg.curveResolution; j++) {
                phi_vals.push(0.1 + (j / (cfg.curveResolution - 1)) * 89.8);
            }
            
            // Use simplified R=1 equation (only cos_phi term)
            const theta_vals = slipLineEquationVecR1(phi_vals, lambda_x, lambda_y, lambda_z, k);
            const valid: boolean[] = theta_vals.map(t => t >= 0.1 && t <= 89.9);
            const validCount = valid.filter(v => v).length;
            console.log(`    Generated ${theta_vals.length} theta values, ${validCount} valid (θ in [0.1, 89.9])°`);
            if (validCount < 5) {
                console.log(`    Skipping: too few valid points`);
                continue;
            }

            const phi_valid: number[] = [];
            const theta_valid: number[] = [];
            for (let i = 0; i < valid.length; i++) {
                if (valid[i]) {
                    phi_valid.push(phi_vals[i]);
                    theta_valid.push(theta_vals[i]);
                }
            }

            const { x, y, z } = sphericalToCartesianVec(phi_valid, theta_valid);
            const normal_stress = calculateNormalStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);

            const reflected = reflectCurveToAllOctants(phi_valid, theta_valid);
            console.log(`    Created ${reflected.length} reflected curves (8 octants)`);
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

        // Add special slip lines: great circles (x,y) and (y,z)
        const specialCurves = generateSpecialSlipLinesR1(lambda_x, lambda_y, lambda_z, cfg.curveResolution);
        console.log(`  Added ${specialCurves.length} special slip lines (great circles)`);
        allCurves.push(...specialCurves);

    // CASE 3: General case (0 < R < 1)
    // Use equidistant points on sigma2 meridian with solved k values
    } else {
        console.log('[Slip Lines] General case');
        
        // Calculate sigma2 meridian azimuth
        const R = (lambda_z - lambda_y) / (lambda_x - lambda_y);
        const phi_sigma2 = Math.acos(Math.sqrt(Math.max(0, Math.min(1, R)))) * 180 / Math.PI;
        
        // Generate equidistant theta points on sigma2 meridian
        const theta_meridian: number[] = [];
        for (let i = 1; i < cfg.nCurves; i++) {
            const t = i / cfg.nCurves;
            theta_meridian.push(t * 90); // θ from 0° to 90°
        }
        
        // Solve for k at each equidistant theta point
        for (const theta_point of theta_meridian) {
            const k = solveForKFromPoint(phi_sigma2, theta_point, lambda_x, lambda_y, lambda_z);
            
            // Generate slip line with this k value, dense sampling near poles
            const phi_vals = [];
            for (let i = 0; i < 50; i++) phi_vals.push(0.5 + (i / 50) * 9.5);
            for (let i = 0; i < cfg.curveResolution - 100; i++) phi_vals.push(10 + (i / (cfg.curveResolution - 100)) * 70);
            for (let i = 0; i < 50; i++) phi_vals.push(80 + (i / 50) * 9.5);
            
            const theta_vals = slipLineEquationVec(phi_vals, lambda_x, lambda_y, lambda_z, k);
            const valid: boolean[] = theta_vals.map(t => t >= 0.1 && t <= 89.9);
            const validCount = valid.filter(v => v).length;
            if (validCount < 5) continue;

            const phi_valid: number[] = [];
            const theta_valid: number[] = [];
            for (let i = 0; i < valid.length; i++) {
                if (valid[i]) {
                    phi_valid.push(phi_vals[i]);
                    theta_valid.push(theta_vals[i]);
                }
            }

            const { x, y, z } = sphericalToCartesianVec(phi_valid, theta_valid);
            const normal_stress = calculateNormalStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_valid, theta_valid, lambda_x, lambda_y, lambda_z);

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

        // Add special slip lines: great circles (x,y), (x,z), and (y,z)
        const specialCurves = generateSpecialSlipLinesGeneral(lambda_x, lambda_y, lambda_z, cfg.curveResolution);
        console.log(`  Added ${specialCurves.length} special slip lines (great circles)`);
        allCurves.push(...specialCurves);
    }

    console.log(`[Slip Lines] Total curves generated: ${allCurves.length}`);
    return allCurves;
}

/**
 * Reflect a special slip line to relevant octants
 * For vertical planes: 4 reflections (specific axes)
 * For horizontal plane: 8 reflections (all octants)
 */
function reflectSpecialSlipLine(
    x: number[],
    y: number[],
    z: number[],
    planeType: 'great_circle_xy' | 'great_circle_xz' | 'great_circle_yz'
): Array<{ x: number[]; y: number[]; z: number[] }> {
    const curves: Array<{ x: number[]; y: number[]; z: number[] }> = [];
    
    if (planeType === 'great_circle_xy') {
        // Horizontal plane (z=0): All 8 octant reflections
        const reflections = [
            { sx: 1, sy: 1, sz: 1 },
            { sx: -1, sy: 1, sz: 1 },
            { sx: 1, sy: -1, sz: 1 },
            { sx: 1, sy: 1, sz: -1 },
            { sx: -1, sy: -1, sz: 1 },
            { sx: -1, sy: 1, sz: -1 },
            { sx: 1, sy: -1, sz: -1 },
            { sx: -1, sy: -1, sz: -1 }
        ];
        
        for (const { sx, sy, sz } of reflections) {
            curves.push({
                x: x.map(v => sx * v),
                y: y.map(v => sy * v),
                z: z.map(v => sz * v)
            });
        }
    } else if (planeType === 'great_circle_xz') {
        // Vertical plane (y=0): 4 reflections along x,z axes only (no y reflection)
        const reflections = [
            { sx: 1, sy: 1, sz: 1 },
            { sx: -1, sy: 1, sz: 1 },
            { sx: 1, sy: 1, sz: -1 },
            { sx: -1, sy: 1, sz: -1 }
        ];
        
        for (const { sx, sy, sz } of reflections) {
            curves.push({
                x: x.map(v => sx * v),
                y: y.map(v => sy * v),
                z: z.map(v => sz * v)
            });
        }
    } else if (planeType === 'great_circle_yz') {
        // Vertical plane (x=0): 4 reflections along y,z axes only (no x reflection)
        const reflections = [
            { sx: 1, sy: 1, sz: 1 },
            { sx: 1, sy: -1, sz: 1 },
            { sx: 1, sy: 1, sz: -1 },
            { sx: 1, sy: -1, sz: -1 }
        ];
        
        for (const { sx, sy, sz } of reflections) {
            curves.push({
                x: x.map(v => sx * v),
                y: y.map(v => sy * v),
                z: z.map(v => sz * v)
            });
        }
    }
    
    return curves;
}

/**
 * Generate special slip lines for general case (0 < R < 1)
 * Three great circles: (x,y), (x,z), and (y,z) planes
 */
function generateSpecialSlipLinesGeneral(
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    resolution: number
): SlipLineCurve[] {
    const curves: SlipLineCurve[] = [];

    // Great circle (x,y) plane: z = 0, φ ∈ [0, 360°], θ = 90°
    const xy_x: number[] = [];
    const xy_y: number[] = [];
    const xy_z: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const phi = (i / resolution) * 360;
        const pt = sphericalToCartesian(phi, 90);
        xy_x.push(pt.x);
        xy_y.push(pt.y);
        xy_z.push(pt.z);
    }

    const xy_normal = calculateNormalStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_shear = calculateShearStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_friction = xy_normal.map((n, i) => Math.abs(n) > 1e-10 ? xy_shear[i] / Math.abs(n) : 0);

    // Reflect (x,y) to 8 octants
    const xy_reflected = reflectSpecialSlipLine(xy_x, xy_y, xy_z, 'great_circle_xy');
    for (const curve of xy_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: xy_normal,
            shear_stress: xy_shear,
            friction: xy_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_xy'
        });
    }

    // Great circle (x,z) plane: y = 0, φ = 0°, varies θ
    const xz_x: number[] = [];
    const xz_y: number[] = [];
    const xz_z: number[] = [];
    const xz_phi: number[] = [];
    const xz_theta: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const theta = (i / (resolution - 1)) * 180;
        const pt = sphericalToCartesian(0, theta);
        xz_x.push(pt.x);
        xz_y.push(pt.y);
        xz_z.push(pt.z);
        xz_phi.push(0);
        xz_theta.push(theta);
    }

    const xz_normal = calculateNormalStressVec(xz_phi, xz_theta, lambda_x, lambda_y, lambda_z);
    const xz_shear = calculateShearStressVec(xz_phi, xz_theta, lambda_x, lambda_y, lambda_z);
    const xz_friction = xz_normal.map((n, i) => Math.abs(n) > 1e-10 ? xz_shear[i] / Math.abs(n) : 0);

    // Reflect (x,z) to 4 octants (x,z reflections only)
    const xz_reflected = reflectSpecialSlipLine(xz_x, xz_y, xz_z, 'great_circle_xz');
    for (const curve of xz_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: xz_normal,
            shear_stress: xz_shear,
            friction: xz_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_xz'
        });
    }

    // Great circle (y,z) plane: x = 0, φ = 90°, varies θ
    const yz_x: number[] = [];
    const yz_y: number[] = [];
    const yz_z: number[] = [];
    const yz_phi: number[] = [];
    const yz_theta: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const theta = (i / (resolution - 1)) * 180;
        const pt = sphericalToCartesian(90, theta);
        yz_x.push(pt.x);
        yz_y.push(pt.y);
        yz_z.push(pt.z);
        yz_phi.push(90);
        yz_theta.push(theta);
    }

    const yz_normal = calculateNormalStressVec(yz_phi, yz_theta, lambda_x, lambda_y, lambda_z);
    const yz_shear = calculateShearStressVec(yz_phi, yz_theta, lambda_x, lambda_y, lambda_z);
    const yz_friction = yz_normal.map((n, i) => Math.abs(n) > 1e-10 ? yz_shear[i] / Math.abs(n) : 0);

    // Reflect (y,z) to 4 octants (y,z reflections only)
    const yz_reflected = reflectSpecialSlipLine(yz_x, yz_y, yz_z, 'great_circle_yz');
    for (const curve of yz_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: yz_normal,
            shear_stress: yz_shear,
            friction: yz_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_yz'
        });
    }

    return curves;
}

/**
 * Generate special slip lines for R=0 case (σ2 = σ3)
 * Two great circles: (x,y) and (x,z) planes
 */
function generateSpecialSlipLinesR0(
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    resolution: number
): SlipLineCurve[] {
    const curves: SlipLineCurve[] = [];

    // Great circle (x,y) plane: z = 0, φ ∈ [0, 360°], θ = 90°
    const xy_x: number[] = [];
    const xy_y: number[] = [];
    const xy_z: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const phi = (i / resolution) * 360;
        const pt = sphericalToCartesian(phi, 90);
        xy_x.push(pt.x);
        xy_y.push(pt.y);
        xy_z.push(pt.z);
    }

    const xy_normal = calculateNormalStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_shear = calculateShearStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_friction = xy_normal.map((n, i) => Math.abs(n) > 1e-10 ? xy_shear[i] / Math.abs(n) : 0);

    // Reflect (x,y) to 8 octants
    const xy_reflected = reflectSpecialSlipLine(xy_x, xy_y, xy_z, 'great_circle_xy');
    for (const curve of xy_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: xy_normal,
            shear_stress: xy_shear,
            friction: xy_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_xy'
        });
    }

    // Great circle (x,z) plane: y = 0, φ = 0°, varies θ
    const xz_x: number[] = [];
    const xz_y: number[] = [];
    const xz_z: number[] = [];
    const xz_phi: number[] = [];
    const xz_theta: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const theta = (i / (resolution - 1)) * 180;
        const pt = sphericalToCartesian(0, theta);
        xz_x.push(pt.x);
        xz_y.push(pt.y);
        xz_z.push(pt.z);
        xz_phi.push(0);
        xz_theta.push(theta);
    }

    const xz_normal = calculateNormalStressVec(xz_phi, xz_theta, lambda_x, lambda_y, lambda_z);
    const xz_shear = calculateShearStressVec(xz_phi, xz_theta, lambda_x, lambda_y, lambda_z);
    const xz_friction = xz_normal.map((n, i) => Math.abs(n) > 1e-10 ? xz_shear[i] / Math.abs(n) : 0);

    // Reflect (x,z) to 4 octants (x,z reflections only)
    const xz_reflected = reflectSpecialSlipLine(xz_x, xz_y, xz_z, 'great_circle_xz');
    for (const curve of xz_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: xz_normal,
            shear_stress: xz_shear,
            friction: xz_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_xz'
        });
    }

    return curves;
}

/**
 * Generate special slip lines for R=1 case (σ2 = σ1)
 * Two great circles: (x,y) and (y,z) planes
 */
function generateSpecialSlipLinesR1(
    lambda_x: number,
    lambda_y: number,
    lambda_z: number,
    resolution: number
): SlipLineCurve[] {
    const curves: SlipLineCurve[] = [];

    // Great circle (x,y) plane: z = 0, φ ∈ [0, 360°], θ = 90°
    const xy_x: number[] = [];
    const xy_y: number[] = [];
    const xy_z: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const phi = (i / resolution) * 360;
        const pt = sphericalToCartesian(phi, 90);
        xy_x.push(pt.x);
        xy_y.push(pt.y);
        xy_z.push(pt.z);
    }

    const xy_normal = calculateNormalStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_shear = calculateShearStressVec(
        Array(resolution).fill(0).map((_, i) => (i / resolution) * 360),
        Array(resolution).fill(90),
        lambda_x, lambda_y, lambda_z
    );
    const xy_friction = xy_normal.map((n, i) => Math.abs(n) > 1e-10 ? xy_shear[i] / Math.abs(n) : 0);

    // Reflect (x,y) to 8 octants
    const xy_reflected = reflectSpecialSlipLine(xy_x, xy_y, xy_z, 'great_circle_xy');
    for (const curve of xy_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: xy_normal,
            shear_stress: xy_shear,
            friction: xy_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_xy'
        });
    }

    // Great circle (y,z) plane: x = 0, φ = 90°, varies θ
    const yz_x: number[] = [];
    const yz_y: number[] = [];
    const yz_z: number[] = [];
    const yz_phi: number[] = [];
    const yz_theta: number[] = [];
    for (let i = 0; i < resolution; i++) {
        const theta = (i / (resolution - 1)) * 180;
        const pt = sphericalToCartesian(90, theta);
        yz_x.push(pt.x);
        yz_y.push(pt.y);
        yz_z.push(pt.z);
        yz_phi.push(90);
        yz_theta.push(theta);
    }

    const yz_normal = calculateNormalStressVec(yz_phi, yz_theta, lambda_x, lambda_y, lambda_z);
    const yz_shear = calculateShearStressVec(yz_phi, yz_theta, lambda_x, lambda_y, lambda_z);
    const yz_friction = yz_normal.map((n, i) => Math.abs(n) > 1e-10 ? yz_shear[i] / Math.abs(n) : 0);

    // Reflect (y,z) to 4 octants (y,z reflections only)
    const yz_reflected = reflectSpecialSlipLine(yz_x, yz_y, yz_z, 'great_circle_yz');
    for (const curve of yz_reflected) {
        curves.push({
            x: curve.x,
            y: curve.y,
            z: curve.z,
            normal_stress: yz_normal,
            shear_stress: yz_shear,
            friction: yz_friction,
            type: 'slip_line',
            isSpecial: true,
            specialType: 'great_circle_yz'
        });
    }

    return curves;
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

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Calculate normal stress at a point on the sphere surface (for coloring)
 * Converts Cartesian to spherical and computes stress value
 */
export function getEquipotentialPoint(
    point: Point3D,
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    const r = Math.sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    if (r < 1e-10) return 0;

    // Convert to spherical
    const theta_rad = Math.acos(Math.max(-1, Math.min(1, point.z / r)));
    let phi_rad = Math.atan2(point.y, point.x);
    if (phi_rad < 0) phi_rad += 2 * Math.PI;

    const theta_deg = (theta_rad * 180) / Math.PI;
    const phi_deg = (phi_rad * 180) / Math.PI;

    return calculateNormalStress(phi_deg, theta_deg, lambda_x, lambda_y, lambda_z);
}