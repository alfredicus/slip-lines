/**
 * Stress Field Computation Module
 * 
 * Implements theoretical slip lines and equipotentials from the Python implementation.
 * Based on: Taboada, A., 1993, Stress and strain from striated pebbles.
 * 
 * IMPROVED SLIP LINE DISTRIBUTION:
 * - Uses sigma2 meridian as reference for even distribution
 * - Correctly handles special great circle slip lines (xy, xz, yz)
 * - Optimized reflections: special lines reflected 4x, general lines reflected 8x
 * - Supports three stress regimes: general (0<R<1), R=0, R=1
 */

// ============================================================================
// TYPES
// ============================================================================

export interface Point3D {
    x: number;
    y: number;
    z: number;
}

/**
 * Identifies which principal plane a special slip line lies on
 * 'xy' = plane where z=0 (sigma1-sigma3 plane)
 * 'xz' = plane where y=0 (sigma1-sigma2 plane)
 * 'yz' = plane where x=0 (sigma2-sigma3 plane)
 */
type SpecialSlipLineType = 'xy' | 'xz' | 'yz' | null;

export interface SlipLineCurve {
    x: number[];
    y: number[];
    z: number[];
    normal_stress: number[];
    shear_stress: number[];
    friction: number[];
    type: 'slip_line' | 'meridian';
    
    // NEW: Metadata for special slip lines (great circles on principal planes)
    isSpecial?: boolean;              // True if this is a great circle
    specialType?: SpecialSlipLineType; // Which plane (xy, xz, yz)
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
// HELPER FUNCTIONS - SPECIAL SLIP LINE IDENTIFICATION
// ============================================================================

/**
 * Determine which special slip lines (great circles) exist for given stress ratio R
 * 
 * CONTEXT:
 * - L-frame axes: L_x=σ1, L_y=σ3, L_z=σ2
 * - R = (σ2 - σ3) / (σ1 - σ3), where R ∈ [0,1]
 * 
 * CASES:
 * - 0 < R < 1 (general): All three great circles exist
 * - R = 0 (σ2=σ3): (xy) and (xz) exist, (yz) is degenerate
 * - R = 1 (σ2=σ1): (xy) and (yz) exist, (xz) is degenerate
 */
function identifySpecialSlipLines(R: number, tolerance: number = 1e-6): SpecialSlipLineType[] {
    const specials: SpecialSlipLineType[] = [];
    
    // (x,y) plane always exists
    specials.push('xy');
    
    if (Math.abs(R) < tolerance) {
        // R = 0: σ2 = σ3
        specials.push('xz');
    } else if (Math.abs(R - 1) < tolerance) {
        // R = 1: σ2 = σ1
        specials.push('yz');
    } else {
        // 0 < R < 1: general case
        specials.push('xz');
        specials.push('yz');
    }
    
    return specials;
}

/**
 * Calculate the azimuth angle phi where the sigma2 meridian lies
 * 
 * The special equipotential is: R = cos²(phi)
 * Therefore: phi = arccos(√R)
 * 
 * @param R - Stress ratio [0,1]
 * @returns Azimuth angle in degrees where sigma2 meridian exists
 */
function calculateSigma2MeridianAzimuth(R: number): number {
    // R = (σ2 - σ3) / (σ1 - σ3) = cos²(phi)
    // Therefore: phi = arccos(√R)
    
    const clampedR = Math.max(0, Math.min(1, R)); // Clamp to [0,1]
    const phi_rad = Math.acos(Math.sqrt(clampedR));
    return (phi_rad * 180) / Math.PI;
}

/**
 * Generate equidistant theta points on the sigma2 meridian
 * 
 * These points will be used to solve for k values, ensuring slip lines
 * are evenly distributed across the sphere surface.
 * 
 * @param nCurves - Number of slip lines to generate
 * @param tolerance - Tolerance for avoiding poles (degrees)
 * @returns Array of theta values in degrees, evenly spaced on meridian
 */
function generateEquidistantPointsOnSigma2Meridian(
    nCurves: number,
    tolerance: number = 0.1
): number[] {
    const thetas: number[] = [];
    const maxTheta = 90;
    
    // Generate nCurves evenly spaced theta values
    for (let i = 1; i < nCurves; i++) {
        const t = i / nCurves; // Parameter from 0 to 1
        const theta = t * maxTheta;
        thetas.push(theta);
    }
    
    return thetas;
}

/**
 * Solve for the integration constant k at a specific point on the slip line
 * 
 * Given a point (phi, theta) that should lie on the slip line, solve for k:
 * 
 * tan(θ) = k × |sin(φ)|^((λ1 - λ2)/(λ3 - λ1)) × |cos(φ)|^((λ2 - λ3)/(λ3 - λ1))
 * 
 * Therefore: k = tan(θ) / [|sin(φ)|^exp_sin × |cos(φ)|^exp_cos]
 * 
 * @param phi_deg - Azimuth angle in degrees
 * @param theta_deg - Colatitude angle in degrees
 * @param lambda_1 - σ1 value (L_x)
 * @param lambda_3 - σ3 value (L_y)
 * @param lambda_2 - σ2 value (L_z)
 * @returns Integration constant k
 */
function solveForKFromPoint(
    phi_deg: number,
    theta_deg: number,
    lambda_1: number,
    lambda_3: number,
    lambda_2: number
): number {
    const phi_rad = (phi_deg * Math.PI) / 180;
    const theta_rad = (theta_deg * Math.PI) / 180;
    
    // Calculate exponents
    const exp_sin = (lambda_1 - lambda_2) / (lambda_3 - lambda_1);
    const exp_cos = (lambda_2 - lambda_3) / (lambda_3 - lambda_1);
    
    // Numerical stability: avoid zero sin/cos
    const sin_phi = Math.max(Math.abs(Math.sin(phi_rad)), 1e-10);
    const cos_phi = Math.max(Math.abs(Math.cos(phi_rad)), 1e-10);
    
    // tan(θ) from the point
    const tan_theta = Math.tan(theta_rad);
    
    // Solve for k
    const denominator = Math.pow(sin_phi, exp_sin) * Math.pow(cos_phi, exp_cos);
    const k = Math.abs(tan_theta) / Math.max(denominator, 1e-10);
    
    return k;
}

/**
 * Generate a single slip line using the parametric equation
 * 
 * Parametric equation: tan(θ) = k × |sin(φ)|^exp_sin × |cos(φ)|^exp_cos
 * 
 * L-frame axes: L_x=σ1, L_y=σ3, L_z=σ2
 * So: λ1=σ1, λ3=σ3, λ2=σ2
 * 
 * @param k - Integration constant
 * @param lambda_1 - σ1 value
 * @param lambda_3 - σ3 value
 * @param lambda_2 - σ2 value
 * @param nPoints - Resolution: number of points on curve
 * @returns Object with phi, theta arrays and Cartesian coordinates
 */
function generateSingleSlipLine(
    k: number,
    lambda_1: number,
    lambda_3: number,
    lambda_2: number,
    nPoints: number = 200
): {
    phi: number[];
    theta: number[];
    x: number[];
    y: number[];
    z: number[];
} {
    const phi_range: number[] = [];
    
    // Dense sampling near poles, sparser in middle
    for (let i = 0; i < 50; i++) {
        phi_range.push(0.5 + (i / 50) * 9.5);
    }
    for (let i = 0; i < nPoints - 100; i++) {
        phi_range.push(10 + (i / (nPoints - 100)) * 70);
    }
    for (let i = 0; i < 50; i++) {
        phi_range.push(80 + (i / 50) * 9.5);
    }
    
    // Generate theta values using slip line equation
    const theta_range: number[] = [];
    for (const phi of phi_range) {
        const theta = slipLineEquation(phi, lambda_1, lambda_3, lambda_2, k);
        if (theta >= 0.1 && theta <= 89.9) {
            theta_range.push(theta);
        }
    }
    
    // Convert to Cartesian
    const { x, y, z } = sphericalToCartesianVec(phi_range, theta_range);
    
    return { phi: phi_range, theta: theta_range, x, y, z };
}

/**
 * Generate a special slip line (great circle on principal plane)
 * 
 * Great circles are generated by setting one coordinate to zero and varying the other two.
 * For example, 'xy' plane has z=0, so we vary x and y.
 * 
 * @param specialType - Which plane: 'xy', 'xz', or 'yz'
 * @param nPoints - Number of points on the curve
 * @returns Object with x, y, z coordinates
 */
function generateSpecialSlipLine(
    specialType: SpecialSlipLineType,
    nPoints: number = 100
): { x: number[]; y: number[]; z: number[] } {
    const x: number[] = [];
    const y: number[] = [];
    const z: number[] = [];
    
    if (specialType === 'xy') {
        // Great circle (x,y) plane: z=0
        // Parametrized as: x = cos(t), y = sin(t), z = 0
        for (let i = 0; i < nPoints; i++) {
            const t = (i / nPoints) * 2 * Math.PI;
            x.push(Math.cos(t));
            y.push(Math.sin(t));
            z.push(0);
        }
    } else if (specialType === 'xz') {
        // Great circle (x,z) plane: y=0
        // Parametrized as: x = cos(t), y = 0, z = sin(t)
        for (let i = 0; i < nPoints; i++) {
            const t = (i / nPoints) * 2 * Math.PI;
            x.push(Math.cos(t));
            y.push(0);
            z.push(Math.sin(t));
        }
    } else if (specialType === 'yz') {
        // Great circle (y,z) plane: x=0
        // Parametrized as: x = 0, y = cos(t), z = sin(t)
        for (let i = 0; i < nPoints; i++) {
            const t = (i / nPoints) * 2 * Math.PI;
            x.push(0);
            y.push(Math.cos(t));
            z.push(Math.sin(t));
        }
    }
    
    return { x, y, z };
}

/**
 * Reflect a great circle to only the relevant octants
 * 
 * OPTIMIZATION: Great circles don't need to be reflected across the axis
 * perpendicular to their plane, because the reflection is identical to the original.
 * 
 * Examples:
 * - (x,y) plane (z=0): Reflect across x,y only (4 reflections, skip z)
 * - (x,z) plane (y=0): Reflect across x,z only (4 reflections, skip y)
 * - (y,z) plane (x=0): Reflect across y,z only (4 reflections, skip x)
 */
function reflectSpecialSlipLineToOctants(
    x: number[],
    y: number[],
    z: number[],
    specialType: SpecialSlipLineType
): Array<{ x: number[]; y: number[]; z: number[] }> {
    const reflections: Array<{ sx: number; sy: number; sz: number }> = [];
    
    if (specialType === 'xy') {
        // z=0: Reflect across x,y only (4 reflections, sz always 1)
        reflections.push({ sx: 1, sy: 1, sz: 1 });
        reflections.push({ sx: -1, sy: 1, sz: 1 });
        reflections.push({ sx: 1, sy: -1, sz: 1 });
        reflections.push({ sx: -1, sy: -1, sz: 1 });
    } else if (specialType === 'xz') {
        // y=0: Reflect across x,z only (4 reflections, sy always 1)
        reflections.push({ sx: 1, sy: 1, sz: 1 });
        reflections.push({ sx: -1, sy: 1, sz: 1 });
        reflections.push({ sx: 1, sy: 1, sz: -1 });
        reflections.push({ sx: -1, sy: 1, sz: -1 });
    } else if (specialType === 'yz') {
        // x=0: Reflect across y,z only (4 reflections, sx always 1)
        reflections.push({ sx: 1, sy: 1, sz: 1 });
        reflections.push({ sx: 1, sy: -1, sz: 1 });
        reflections.push({ sx: 1, sy: 1, sz: -1 });
        reflections.push({ sx: 1, sy: -1, sz: -1 });
    }
    
    const curves: Array<{ x: number[]; y: number[]; z: number[] }> = [];
    
    for (const { sx, sy, sz } of reflections) {
        const x_refl = x.map((v) => sx * v);
        const y_refl = y.map((v) => sy * v);
        const z_refl = z.map((v) => sz * v);
        curves.push({ x: x_refl, y: y_refl, z: z_refl });
    }
    
    return curves;
}

// ============================================================================
// CURVE REFLECTION (GENERAL - 8 OCTANTS)
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
// SLIP LINE GENERATION - MAIN ALGORITHM
// ============================================================================

/**
 * Generate slip lines with improved distribution using sigma2 meridian reference
 * 
 * ALGORITHM:
 * 1. Calculate stress ratio R from sigma values
 * 2. Identify which special slip lines (great circles) exist for this stress regime
 * 3. Calculate sigma2 meridian azimuth angle
 * 4. Generate evenly-spaced points on the sigma2 meridian
 * 5. At each point, solve for the integration constant k
 * 6. Generate slip line using this k value
 * 7. Reflect each slip line appropriately (4x for specials, 8x for general)
 * 
 * STRESS REGIMES:
 * - General (0 < R < 1): All three special slip lines included
 * - R = 0 (σ2=σ3): Special lines (xy), (xz); revolution meridians
 * - R = 1 (σ2=σ1): Special lines (xy), (yz); revolution meridians
 */
export function generateSlipLines(
    lambda_x: number,   // L_x = σ1 (minimum/most compressive)
    lambda_y: number,   // L_y = σ3 (maximum/most extensional)
    lambda_z: number,   // L_z = σ2 (intermediate)
    config?: Partial<ComputeConfig>
): SlipLineCurve[] {
    const cfg = {
        nCurves: config?.nCurves ?? DEFAULT_N_CURVES,
        curveResolution: config?.curveResolution ?? DEFAULT_CURVE_RESOLUTION
    };

    const allCurves: SlipLineCurve[] = [];

    // ========================================================================
    // STEP 1: CALCULATE STRESS RATIO R
    // ========================================================================
    // R = (σ2 - σ3) / (σ1 - σ3) = (lambda_z - lambda_y) / (lambda_x - lambda_y)
    
    const denominator = lambda_x - lambda_y;
    const R = Math.abs(denominator) > 1e-10 ? (lambda_z - lambda_y) / denominator : 0.5;
    
    // ========================================================================
    // CASE 1: GENERAL CASE (0 < R < 1)
    // ========================================================================
    
    if (Math.abs(R) > 1e-6 && Math.abs(R - 1) > 1e-6) {
        console.log(`[Slip Lines] General case: R = ${R.toFixed(3)}`);
        
        // STEP 2: Calculate sigma2 meridian azimuth
        const phi_sigma2 = calculateSigma2MeridianAzimuth(R);
        console.log(`[Slip Lines] Sigma2 meridian at phi = ${phi_sigma2.toFixed(2)}°`);
        
        // STEP 3: Generate evenly-spaced theta points on sigma2 meridian
        const theta_equidistant = generateEquidistantPointsOnSigma2Meridian(cfg.nCurves);
        console.log(`[Slip Lines] Generated ${cfg.nCurves} equidistant points on sigma2 meridian`);
        
        // STEP 4: Solve for k values at each equidistant point
        const k_values: number[] = [];
        for (const theta of theta_equidistant) {
            const k = solveForKFromPoint(phi_sigma2, theta, lambda_x, lambda_y, lambda_z);
            k_values.push(k);
        }
        
        // STEP 5: Generate slip lines for each k value
        for (const k of k_values) {
            const slipLine = generateSingleSlipLine(k, lambda_x, lambda_y, lambda_z, cfg.curveResolution);
            
            // Calculate stress properties
            const normal_stress = calculateNormalStressVec(
                slipLine.phi, slipLine.theta, lambda_x, lambda_y, lambda_z
            );
            const shear_stress = calculateShearStressVec(
                slipLine.phi, slipLine.theta, lambda_x, lambda_y, lambda_z
            );
            const friction = calculateFrictionVec(
                slipLine.phi, slipLine.theta, lambda_x, lambda_y, lambda_z
            );
            
            // STEP 6: Reflect to all 8 octants
            const reflected = reflectCurveToAllOctants(slipLine.phi, slipLine.theta);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'slip_line',
                    isSpecial: false
                });
            }
        }
        
        // STEP 7: Add special great circle slip lines
        const specials = identifySpecialSlipLines(R);
        console.log(`[Slip Lines] Adding special slip lines: ${specials.join(', ')}`);
        
        for (const specialType of specials) {
            const special = generateSpecialSlipLine(specialType, cfg.curveResolution);
            
            // Calculate stress properties by converting to spherical coordinates
            const phi_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.atan2(special.y[i], special.x[i]) * 180 / Math.PI
            );
            const theta_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.acos(Math.max(-1, Math.min(1, special.z[i]))) * 180 / Math.PI
            );
            
            const normal_stress = calculateNormalStressVec(
                phi_special, theta_special, lambda_x, lambda_y, lambda_z
            );
            const shear_stress = calculateShearStressVec(
                phi_special, theta_special, lambda_x, lambda_y, lambda_z
            );
            const friction = calculateFrictionVec(
                phi_special, theta_special, lambda_x, lambda_y, lambda_z
            );
            
            // Reflect only to relevant octants (4 reflections, not 8)
            const reflected = reflectSpecialSlipLineToOctants(special.x, special.y, special.z, specialType);
            console.log(`[Slip Lines] Special '${specialType}': 4 reflections applied (50% optimization)`);
            
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'slip_line',
                    isSpecial: true,
                    specialType
                });
            }
        }
    }
    
    // ========================================================================
    // CASE 2: REVOLUTION CASE R = 0 (σ2 = σ3)
    // ========================================================================
    
    else if (Math.abs(R) < 1e-6) {
        console.log(`[Slip Lines] Revolution case: R = 0 (σ2 = σ3)`);
        
        // Meridians evenly spaced around (y,z) plane (phi = 90°)
        for (let i = 0; i < cfg.nCurves; i++) {
            const phi = 90; // Fixed meridian at phi = 90
            const theta_vals: number[] = [];
            for (let j = 0; j < cfg.curveResolution; j++) {
                theta_vals.push(0.5 + (j / (cfg.curveResolution - 1)) * 179);
            }
            
            const phi_vals = theta_vals.map(() => phi);
            const { x, y, z } = sphericalToCartesianVec(phi_vals, theta_vals);
            const normal_stress = calculateNormalStressVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            
            const reflected = reflectCurveToAllOctants(phi_vals, theta_vals);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'meridian',
                    isSpecial: false
                });
            }
        }
        
        // Add special slip lines (xy) and (xz), excluding (yz)
        console.log(`[Slip Lines] Adding special slip lines: xy, xz (excluding yz)`);
        
        for (const specialType of ['xy', 'xz'] as SpecialSlipLineType[]) {
            const special = generateSpecialSlipLine(specialType, cfg.curveResolution);
            
            const phi_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.atan2(special.y[i], special.x[i]) * 180 / Math.PI
            );
            const theta_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.acos(Math.max(-1, Math.min(1, special.z[i]))) * 180 / Math.PI
            );
            
            const normal_stress = calculateNormalStressVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            
            const reflected = reflectSpecialSlipLineToOctants(special.x, special.y, special.z, specialType);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'slip_line',
                    isSpecial: true,
                    specialType
                });
            }
        }
    }
    
    // ========================================================================
    // CASE 3: REVOLUTION CASE R = 1 (σ2 = σ1)
    // ========================================================================
    
    else if (Math.abs(R - 1) < 1e-6) {
        console.log(`[Slip Lines] Revolution case: R = 1 (σ2 = σ1)`);
        
        // Meridians evenly spaced around (x,z) plane (phi = 0°)
        for (let i = 0; i < cfg.nCurves; i++) {
            const phi = 0; // Fixed meridian at phi = 0
            const theta_vals: number[] = [];
            for (let j = 0; j < cfg.curveResolution; j++) {
                theta_vals.push(0.5 + (j / (cfg.curveResolution - 1)) * 179);
            }
            
            const phi_vals = theta_vals.map(() => phi);
            const { x, y, z } = sphericalToCartesianVec(phi_vals, theta_vals);
            const normal_stress = calculateNormalStressVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_vals, theta_vals, lambda_x, lambda_y, lambda_z);
            
            const reflected = reflectCurveToAllOctants(phi_vals, theta_vals);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'meridian',
                    isSpecial: false
                });
            }
        }
        
        // Add special slip lines (xy) and (yz), excluding (xz)
        console.log(`[Slip Lines] Adding special slip lines: xy, yz (excluding xz)`);
        
        for (const specialType of ['xy', 'yz'] as SpecialSlipLineType[]) {
            const special = generateSpecialSlipLine(specialType, cfg.curveResolution);
            
            const phi_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.atan2(special.y[i], special.x[i]) * 180 / Math.PI
            );
            const theta_special = Array(special.x.length).fill(0).map((_, i) => 
                Math.acos(Math.max(-1, Math.min(1, special.z[i]))) * 180 / Math.PI
            );
            
            const normal_stress = calculateNormalStressVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            const shear_stress = calculateShearStressVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            const friction = calculateFrictionVec(phi_special, theta_special, lambda_x, lambda_y, lambda_z);
            
            const reflected = reflectSpecialSlipLineToOctants(special.x, special.y, special.z, specialType);
            for (const curve of reflected) {
                allCurves.push({
                    x: curve.x,
                    y: curve.y,
                    z: curve.z,
                    normal_stress,
                    shear_stress,
                    friction,
                    type: 'slip_line',
                    isSpecial: true,
                    specialType
                });
            }
        }
    }

    console.log(`[Slip Lines] Total curves generated: ${allCurves.length}`);
    return allCurves;
}

/**
 * Get equipotential value at a point on sphere
 */
export function getEquipotentialPoint(
    point: { x: number; y: number; z: number },
    lambda_x: number,
    lambda_y: number,
    lambda_z: number
): number {
    return (
        (lambda_x * point.x * point.x + lambda_y * point.y * point.y) +
        lambda_z * point.z * point.z
    );
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