import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls';
import {
    generateEquipotentials,
    generateSlipLines,
    generatePrincipalAxes,
    EquipotentialCurve,
    SlipLineCurve,
    PrincipalAxis,
} from './stressFieldCompute';

// ============================================================================
// TYPE DEFINITIONS
// ============================================================================

interface ControlState {
    sigmaX: number;
    sigmaY: number;
    sigmaZ: number;
    showEquipotential: boolean;
    showSlipLines: boolean;
    showStressAxes: boolean;
    slipLineColor: string;
    equipotentialColor: string;
}

// ============================================================================
// CONSTANTS
// ============================================================================

const INITIAL_CONTROLS: ControlState = {
    sigmaX: 1.0,
    sigmaY: 3.0,
    sigmaZ: 2.0,
    showEquipotential: true,
    showSlipLines: true,
    showStressAxes: true,
    slipLineColor: '#D62728',
    equipotentialColor: '#1F77B4'
};

const COMPUTE_CONFIG = {
    nCurves: 12,
    nIsoCurves: 8,
    curveResolution: 200
};

const CAMERA_INITIAL_POS = { x: 3, y: 3, z: 3 };
const SPHERE_OPACITY = 1;
const AXES_SIZE = 0.8;

// ============================================================================
// REACT COMPONENT
// ============================================================================

const SlipLinesVisualization: React.FC = () => {
    // -------------------------------------------------------------------------
    // REFS
    // -------------------------------------------------------------------------

    const containerRef = useRef<HTMLDivElement>(null);
    const sceneRef = useRef<THREE.Scene | null>(null);
    const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
    const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
    const sphereGroupRef = useRef<THREE.Group>(new THREE.Group());

    // -------------------------------------------------------------------------
    // STATE
    // -------------------------------------------------------------------------

    const [controls, setControls] = useState<ControlState>(INITIAL_CONTROLS);

    // -------------------------------------------------------------------------
    // EVENT HANDLERS
    // -------------------------------------------------------------------------

    const handleControlChange = <K extends keyof ControlState>(
        key: K,
        value: ControlState[K]
    ): void => {
        setControls((prev) => ({ ...prev, [key]: value }));
    };

    // -------------------------------------------------------------------------
    // THREE.JS SCENE INITIALIZATION
    // -------------------------------------------------------------------------

    useEffect(() => {
        if (!containerRef.current) return;

        // Create scene
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0xf0f0f0);
        sceneRef.current = scene;

        // Create camera
        const width = containerRef.current.clientWidth;
        const height = containerRef.current.clientHeight;
        const camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        camera.position.set(CAMERA_INITIAL_POS.x, CAMERA_INITIAL_POS.y, CAMERA_INITIAL_POS.z);
        camera.lookAt(0, 0, 0);
        cameraRef.current = camera;

        // Create renderer
        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(width, height);
        renderer.setPixelRatio(window.devicePixelRatio);
        containerRef.current.appendChild(renderer.domElement);
        rendererRef.current = renderer;

        // -----------------------------------------------------------------------
        // LIGHTING
        // -----------------------------------------------------------------------

        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        scene.add(ambientLight);

        const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight1.position.set(5, 5, 5);
        scene.add(directionalLight1);

        const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight2.position.set(-5, -5, -5);
        scene.add(directionalLight2);

        // -----------------------------------------------------------------------
        // SPHERE MESH
        // -----------------------------------------------------------------------

        const sphereGeometry = new THREE.SphereGeometry(1, 64, 64);
        const sphereMaterial = new THREE.MeshPhongMaterial({
            color: 0xe8e8e8,
            wireframe: false,
            opacity: SPHERE_OPACITY,
            transparent: true
        });
        const sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
        scene.add(sphere);

        // -----------------------------------------------------------------------
        // COORDINATE AXES HELPER
        // -----------------------------------------------------------------------

        const axesHelper = new THREE.AxesHelper(AXES_SIZE);
        scene.add(axesHelper);

        // Add the group for curves
        scene.add(sphereGroupRef.current);

        // -----------------------------------------------------------------------
        // TRACKBALL CONTROLS
        // -----------------------------------------------------------------------

        const controls = new TrackballControls(camera, renderer.domElement);
        controls.rotateSpeed = 1.0;
        controls.zoomSpeed = 1.2;
        controls.panSpeed = 0.8;
        controls.noRotate = false;
        controls.noZoom = false;
        controls.noPan = false;
        controls.autoRotate = false;
        controls.target.set(0, 0, 0);

        // -----------------------------------------------------------------------
        // ANIMATION LOOP
        // -----------------------------------------------------------------------

        let animationId: number;
        const animate = (): void => {
            animationId = requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        };
        animate();

        // -----------------------------------------------------------------------
        // WINDOW RESIZE HANDLER
        // -----------------------------------------------------------------------

        const handleResize = (): void => {
            if (!containerRef.current || !cameraRef.current) return;

            const newWidth = containerRef.current.clientWidth;
            const newHeight = containerRef.current.clientHeight;

            cameraRef.current.aspect = newWidth / newHeight;
            cameraRef.current.updateProjectionMatrix();
            renderer.setSize(newWidth, newHeight);
        };

        window.addEventListener('resize', handleResize);

        // -----------------------------------------------------------------------
        // CLEANUP
        // -----------------------------------------------------------------------

        return () => {
            window.removeEventListener('resize', handleResize);
            
            controls.dispose();
            cancelAnimationFrame(animationId);
            renderer.dispose();

            if (containerRef.current?.contains(renderer.domElement)) {
                containerRef.current.removeChild(renderer.domElement);
            }
        };
    }, []);

    // =========================================================================
    // VISUALIZATION UPDATE - Compute and render curves
    // =========================================================================

    useEffect(() => {
        if (!sphereGroupRef.current) return;

        // Clear previous geometry
        sphereGroupRef.current.clear();

        // -----------------------------------------------------------------------
        // GENERATE AND ADD EQUIPOTENTIALS
        // -----------------------------------------------------------------------

        if (controls.showEquipotential) {
            const equipotentials = generateEquipotentials(
                controls.sigmaX,
                controls.sigmaY,
                controls.sigmaZ,
                COMPUTE_CONFIG
            );
            addLinesToScene(equipotentials, controls.equipotentialColor, sphereGroupRef.current);
        }

        // -----------------------------------------------------------------------
        // GENERATE AND ADD SLIP LINES
        // -----------------------------------------------------------------------

        if (controls.showSlipLines) {
            const slipLines = generateSlipLines(
                controls.sigmaX,
                controls.sigmaY,
                controls.sigmaZ,
                COMPUTE_CONFIG
            );
            addLinesToScene(slipLines, controls.slipLineColor, sphereGroupRef.current);
        }

        // -----------------------------------------------------------------------
        // GENERATE AND ADD PRINCIPAL STRESS AXES
        // -----------------------------------------------------------------------

        if (controls.showStressAxes) {
            const axes = generatePrincipalAxes(
                controls.sigmaX,
                controls.sigmaY,
                controls.sigmaZ
            );
            addAxesToScene(axes, sphereGroupRef.current);
        }
    }, [
        controls.sigmaX,
        controls.sigmaY,
        controls.sigmaZ,
        controls.showEquipotential,
        controls.showSlipLines,
        controls.showStressAxes,
        controls.slipLineColor,
        controls.equipotentialColor
    ]);

    // =========================================================================
    // HELPER FUNCTIONS - THREE.JS GEOMETRY CREATION
    // =========================================================================

    /**
     * Add curve lines to the Three.js scene
     * Handles both equipotential curves and slip lines
     */
    const addLinesToScene = (
        curves: EquipotentialCurve[] | SlipLineCurve[],
        color: string,
        group: THREE.Group
    ): void => {
        (curves as any[]).forEach((curve) => {
            // Extract coordinates from curve
            const x = (curve as any).x || [];
            const y = (curve as any).y || [];
            const z = (curve as any).z || [];

            if (x.length === 0 || y.length === 0 || z.length === 0) return;
            if (x.length !== y.length || y.length !== z.length) return;

            // Create Three.js geometry
            const geometry = new THREE.BufferGeometry();
            const positions: number[] = [];

            for (let i = 0; i < x.length; i++) {
                positions.push(x[i], y[i], z[i]);
            }

            geometry.setAttribute(
                'position',
                new THREE.BufferAttribute(new Float32Array(positions), 3)
            );

            // Create material and line
            const material = new THREE.LineBasicMaterial({
                color: new THREE.Color(color),
                fog: false
            });

            const lineObj = new THREE.Line(geometry, material);
            group.add(lineObj);
        });
    };

    /**
     * Add principal stress axes to the scene with arrow heads
     */
    const addAxesToScene = (axes: PrincipalAxis[], group: THREE.Group): void => {
        axes.forEach((axis) => {
            // ---------------------------------------------------------------
            // MAIN AXIS LINE
            // ---------------------------------------------------------------

            const geometry = new THREE.BufferGeometry();
            const positions = [
                0, 0, 0,
                axis.axis.x * 1.5, axis.axis.y * 1.5, axis.axis.z * 1.5
            ];

            geometry.setAttribute(
                'position',
                new THREE.BufferAttribute(new Float32Array(positions), 3)
            );

            const material = new THREE.LineBasicMaterial({ color: axis.colorHex });
            const line = new THREE.Line(geometry, material);
            group.add(line);

            // ---------------------------------------------------------------
            // ARROW HEAD
            // ---------------------------------------------------------------

            const arrowDir = new THREE.Vector3(axis.axis.x, axis.axis.y, axis.axis.z).normalize();
            const arrowLength = 0.2;
            const arrowOrigin = new THREE.Vector3(
                axis.axis.x * 1.5,
                axis.axis.y * 1.5,
                axis.axis.z * 1.5
            );

            // Get perpendicular vectors for arrow wings
            let perpDir1 = new THREE.Vector3()
                .crossVectors(arrowDir, new THREE.Vector3(0, 0, 1))
                .normalize();

            if (perpDir1.length() === 0) {
                perpDir1.set(1, 0, 0);
            }

            const perpDir2 = new THREE.Vector3()
                .crossVectors(arrowDir, perpDir1)
                .normalize();

            // Create arrow head geometry
            const arrowGeometry = new THREE.BufferGeometry();
            const arrowPositions = [
                // First wing
                arrowOrigin.x, arrowOrigin.y, arrowOrigin.z,
                arrowOrigin.x - arrowDir.x * arrowLength + perpDir1.x * arrowLength * 0.3,
                arrowOrigin.y - arrowDir.y * arrowLength + perpDir1.y * arrowLength * 0.3,
                arrowOrigin.z - arrowDir.z * arrowLength + perpDir1.z * arrowLength * 0.3,
                // Second wing
                arrowOrigin.x, arrowOrigin.y, arrowOrigin.z,
                arrowOrigin.x - arrowDir.x * arrowLength + perpDir2.x * arrowLength * 0.3,
                arrowOrigin.y - arrowDir.y * arrowLength + perpDir2.y * arrowLength * 0.3,
                arrowOrigin.z - arrowDir.z * arrowLength + perpDir2.z * arrowLength * 0.3
            ];

            arrowGeometry.setAttribute(
                'position',
                new THREE.BufferAttribute(new Float32Array(arrowPositions), 3)
            );

            const arrow = new THREE.LineSegments(arrowGeometry, material);
            group.add(arrow);
        });
    };

    // =========================================================================
    // RENDER
    // =========================================================================

    return (
        <div style={{ display: 'flex', height: '100vh', fontFamily: 'Arial, sans-serif' }}>
            {/* 3D Visualization */}
            <div
                ref={containerRef}
                style={{
                    flex: 1,
                    position: 'relative'
                }}
            />

            {/* Control Panel */}
            <ControlPanel controls={controls} onControlChange={handleControlChange} />
        </div>
    );
};

// ============================================================================
// CONTROL PANEL COMPONENT
// ============================================================================

interface ControlPanelProps {
    controls: ControlState;
    onControlChange: <K extends keyof ControlState>(key: K, value: ControlState[K]) => void;
}

const ControlPanel: React.FC<ControlPanelProps> = ({ controls, onControlChange }) => {
    return (
        <div
            style={{
                width: 350,
                padding: 20,
                backgroundColor: '#f5f5f5',
                borderLeft: '1px solid #ddd',
                overflowY: 'auto'
            }}
        >
            <h2 style={{ marginTop: 0, fontSize: 18, fontWeight: 'bold' }}>
                3D Slip Lines & Equipotentials
            </h2>

            {/* Stress Tensor Parameters Section */}
            <div style={{ marginBottom: 25 }}>
                <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
                    Stress Tensor (σ)
                </h3>

                <StressSlider
                    label="σ_x (σ₁ - minimum)"
                    value={controls.sigmaX}
                    onChange={(val) => onControlChange('sigmaX', val)}
                />

                <StressSlider
                    label="σ_y (σ₃ - maximum)"
                    value={controls.sigmaY}
                    onChange={(val) => onControlChange('sigmaY', val)}
                />

                <StressSlider
                    label="σ_z (σ₂ - intermediate)"
                    value={controls.sigmaZ}
                    onChange={(val) => onControlChange('sigmaZ', val)}
                />
            </div>

            {/* Visualization Options Section */}
            <div style={{ marginBottom: 25 }}>
                <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
                    Visualization
                </h3>

                <ToggleControl
                    label="Show Equipotentials"
                    checked={controls.showEquipotential}
                    onChange={(val) => onControlChange('showEquipotential', val)}
                >
                    <ColorPicker
                        value={controls.equipotentialColor}
                        onChange={(val) => onControlChange('equipotentialColor', val)}
                    />
                </ToggleControl>

                <ToggleControl
                    label="Show Slip Lines"
                    checked={controls.showSlipLines}
                    onChange={(val) => onControlChange('showSlipLines', val)}
                >
                    <ColorPicker
                        value={controls.slipLineColor}
                        onChange={(val) => onControlChange('slipLineColor', val)}
                    />
                </ToggleControl>

                <ToggleControl
                    label="Show Stress Axes"
                    checked={controls.showStressAxes}
                    onChange={(val) => onControlChange('showStressAxes', val)}
                />
            </div>

            {/* Info Section */}
            <div style={{
                fontSize: 12,
                color: '#666',
                borderTop: '1px solid #ddd',
                paddingTop: 15
            }}>
                <p style={{ margin: '0 0 8px 0', fontWeight: 'bold' }}>Controls:</p>
                <ul style={{ margin: 0, paddingLeft: 16, lineHeight: 1.6 }}>
                    <li>Left drag: Rotate</li>
                    <li>Middle drag: Zoom</li>
                    <li>Right drag: Pan</li>
                    <li>Scroll: Zoom</li>
                </ul>
            </div>
        </div>
    );
};

// ============================================================================
// UI CONTROL COMPONENTS
// ============================================================================

interface StressSliderProps {
    label: string;
    value: number;
    onChange: (value: number) => void;
}

const StressSlider: React.FC<StressSliderProps> = ({ label, value, onChange }) => (
    <div style={{ marginBottom: 15 }}>
        <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
            {label}: {value.toFixed(2)}
        </label>
        <input
            type="range"
            min="-5"
            max="5"
            step="0.1"
            value={value}
            onChange={(e) => onChange(parseFloat(e.target.value))}
            style={{ width: '100%', cursor: 'pointer' }}
        />
    </div>
);

interface ToggleControlProps {
    label: string;
    checked: boolean;
    onChange: (checked: boolean) => void;
    children?: React.ReactNode;
}

const ToggleControl: React.FC<ToggleControlProps> = ({ label, checked, onChange, children }) => (
    <>
        <div style={{ marginBottom: 12 }}>
            <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer', fontSize: 13 }}>
                <input
                    type="checkbox"
                    checked={checked}
                    onChange={(e) => onChange(e.target.checked)}
                    style={{ marginRight: 8, width: 16, height: 16, cursor: 'pointer' }}
                />
                <span>{label}</span>
            </label>
        </div>
        {checked && children && <div style={{ marginBottom: 12, paddingLeft: 24 }}>{children}</div>}
    </>
);

interface ColorPickerProps {
    value: string;
    onChange: (value: string) => void;
}

const ColorPicker: React.FC<ColorPickerProps> = ({ value, onChange }) => (
    <label style={{ display: 'flex', alignItems: 'center', gap: 8, fontSize: 13 }}>
        Color:
        <input
            type="color"
            value={value}
            onChange={(e) => onChange(e.target.value)}
            style={{ width: 40, height: 30, cursor: 'pointer', border: 'none' }}
        />
    </label>
);

export default SlipLinesVisualization;
// import React, { useEffect, useRef, useState } from 'react';
// import * as THREE from 'three';
// import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls';
// import {
//     generateEquipotentials,
//     generateSlipLines,
//     generatePrincipalAxes,
//     EquipotentialCurve,
//     SlipLineCurve,
//     PrincipalAxis,
// } from './stressFieldCompute';

// // ============================================================================
// // TYPE DEFINITIONS
// // ============================================================================

// interface ControlState {
//     sigmaX: number;
//     sigmaY: number;
//     sigmaZ: number;
//     showEquipotential: boolean;
//     showSlipLines: boolean;
//     showStressAxes: boolean;
//     slipLineColor: string;
//     equipotentialColor: string;
// }

// // ============================================================================
// // CONSTANTS
// // ============================================================================

// const INITIAL_CONTROLS: ControlState = {
//     sigmaX: 1.0,
//     sigmaY: 3.0,
//     sigmaZ: 2.0,
//     showEquipotential: true,
//     showSlipLines: true,
//     showStressAxes: true,
//     slipLineColor: '#D62728',
//     equipotentialColor: '#1F77B4'
// };

// const COMPUTE_CONFIG = {
//     nCurves: 12,
//     nIsoCurves: 8,
//     curveResolution: 200
// };

// const CAMERA_INITIAL_POS = { x: 3, y: 3, z: 3 };
// const SPHERE_OPACITY = 1;
// const AXES_SIZE = 0.8;

// // ============================================================================
// // REACT COMPONENT
// // ============================================================================

// const SlipLinesVisualization: React.FC = () => {
//     // -------------------------------------------------------------------------
//     // REFS
//     // -------------------------------------------------------------------------

//     const containerRef = useRef<HTMLDivElement>(null);
//     const sceneRef = useRef<THREE.Scene | null>(null);
//     const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
//     const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
//     const sphereGroupRef = useRef<THREE.Group>(new THREE.Group());

//     // -------------------------------------------------------------------------
//     // STATE
//     // -------------------------------------------------------------------------

//     const [controls, setControls] = useState<ControlState>(INITIAL_CONTROLS);

//     // -------------------------------------------------------------------------
//     // EVENT HANDLERS
//     // -------------------------------------------------------------------------

//     const handleControlChange = <K extends keyof ControlState>(
//         key: K,
//         value: ControlState[K]
//     ): void => {
//         setControls((prev) => ({ ...prev, [key]: value }));
//     };

//     // -------------------------------------------------------------------------
//     // THREE.JS SCENE INITIALIZATION
//     // -------------------------------------------------------------------------

//     useEffect(() => {
//         if (!containerRef.current) return;

//         // Create scene
//         const scene = new THREE.Scene();
//         scene.background = new THREE.Color(0xf0f0f0);
//         sceneRef.current = scene;

//         // Create camera
//         const width = containerRef.current.clientWidth;
//         const height = containerRef.current.clientHeight;
//         const camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
//         camera.position.set(CAMERA_INITIAL_POS.x, CAMERA_INITIAL_POS.y, CAMERA_INITIAL_POS.z);
//         camera.lookAt(0, 0, 0);
//         cameraRef.current = camera;

//         // Create renderer
//         const renderer = new THREE.WebGLRenderer({ antialias: true });
//         renderer.setSize(width, height);
//         renderer.setPixelRatio(window.devicePixelRatio);
//         containerRef.current.appendChild(renderer.domElement);
//         rendererRef.current = renderer;

//         // -----------------------------------------------------------------------
//         // LIGHTING
//         // -----------------------------------------------------------------------

//         const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
//         scene.add(ambientLight);

//         const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
//         directionalLight.position.set(5, 5, 5);
//         scene.add(directionalLight);

//         // -----------------------------------------------------------------------
//         // SPHERE MESH
//         // -----------------------------------------------------------------------

//         const sphereGeometry = new THREE.SphereGeometry(1, 64, 64);
//         const sphereMaterial = new THREE.MeshPhongMaterial({
//             color: 0xe8e8e8,
//             wireframe: false,
//             opacity: SPHERE_OPACITY,
//             transparent: true
//         });
//         const sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
//         scene.add(sphere);

//         // -----------------------------------------------------------------------
//         // COORDINATE AXES HELPER
//         // -----------------------------------------------------------------------

//         const axesHelper = new THREE.AxesHelper(AXES_SIZE);
//         scene.add(axesHelper);

//         // Add the group for curves
//         scene.add(sphereGroupRef.current);

//         // -----------------------------------------------------------------------
//         // MOUSE CONTROLS
//         // -----------------------------------------------------------------------

//         const mouseState = { down: false, x: 0, y: 0 };

//         const onMouseDown = (e: MouseEvent): void => {
//             mouseState.down = true;
//             mouseState.x = e.clientX;
//             mouseState.y = e.clientY;
//         };

//         const onMouseMove = (e: MouseEvent): void => {
//             if (!mouseState.down || !cameraRef.current) return;

//             const deltaX = e.clientX - mouseState.x;
//             const deltaY = e.clientY - mouseState.y;

//             // Rotate camera around origin
//             cameraRef.current.position.applyAxisAngle(new THREE.Vector3(0, 1, 0), deltaX * 0.01);
//             cameraRef.current.position.applyAxisAngle(new THREE.Vector3(1, 0, 0), deltaY * 0.01);
//             cameraRef.current.lookAt(0, 0, 0);

//             mouseState.x = e.clientX;
//             mouseState.y = e.clientY;
//         };

//         const onMouseUp = (): void => {
//             mouseState.down = false;
//         };

//         const onWheel = (e: WheelEvent): void => {
//             if (!cameraRef.current) return;
//             e.preventDefault();

//             const zoom = e.deltaY > 0 ? 1.1 : 0.9;
//             cameraRef.current.position.multiplyScalar(zoom);
//         };

//         renderer.domElement.addEventListener('mousedown', onMouseDown);
//         renderer.domElement.addEventListener('mousemove', onMouseMove);
//         renderer.domElement.addEventListener('mouseup', onMouseUp);
//         renderer.domElement.addEventListener('wheel', onWheel, { passive: false });

//         // -----------------------------------------------------------------------
//         // ANIMATION LOOP
//         // -----------------------------------------------------------------------

//         let animationId: number;
//         const animate = (): void => {
//             animationId = requestAnimationFrame(animate);
//             renderer.render(scene, camera);
//         };
//         animate();

//         // -----------------------------------------------------------------------
//         // WINDOW RESIZE HANDLER
//         // -----------------------------------------------------------------------

//         const handleResize = (): void => {
//             if (!containerRef.current || !cameraRef.current) return;

//             const newWidth = containerRef.current.clientWidth;
//             const newHeight = containerRef.current.clientHeight;

//             cameraRef.current.aspect = newWidth / newHeight;
//             cameraRef.current.updateProjectionMatrix();
//             renderer.setSize(newWidth, newHeight);
//         };

//         window.addEventListener('resize', handleResize);

//         // -----------------------------------------------------------------------
//         // CLEANUP
//         // -----------------------------------------------------------------------

//         return () => {
//             window.removeEventListener('resize', handleResize);
//             renderer.domElement.removeEventListener('mousedown', onMouseDown);
//             renderer.domElement.removeEventListener('mousemove', onMouseMove);
//             renderer.domElement.removeEventListener('mouseup', onMouseUp);
//             renderer.domElement.removeEventListener('wheel', onWheel);

//             cancelAnimationFrame(animationId);
//             renderer.dispose();

//             if (containerRef.current?.contains(renderer.domElement)) {
//                 containerRef.current.removeChild(renderer.domElement);
//             }
//         };
//     }, []);

//     // =========================================================================
//     // VISUALIZATION UPDATE - Compute and render curves
//     // =========================================================================

//     useEffect(() => {
//         if (!sphereGroupRef.current) return;

//         // Clear previous geometry
//         sphereGroupRef.current.clear();

//         // -----------------------------------------------------------------------
//         // GENERATE AND ADD EQUIPOTENTIALS
//         // -----------------------------------------------------------------------

//         if (controls.showEquipotential) {
//             const equipotentials = generateEquipotentials(
//                 controls.sigmaX,
//                 controls.sigmaY,
//                 controls.sigmaZ,
//                 COMPUTE_CONFIG
//             );
//             addLinesToScene(equipotentials, controls.equipotentialColor, sphereGroupRef.current);
//         }

//         // -----------------------------------------------------------------------
//         // GENERATE AND ADD SLIP LINES
//         // -----------------------------------------------------------------------

//         if (controls.showSlipLines) {
//             const slipLines = generateSlipLines(
//                 controls.sigmaX,
//                 controls.sigmaY,
//                 controls.sigmaZ,
//                 COMPUTE_CONFIG
//             );
//             addLinesToScene(slipLines, controls.slipLineColor, sphereGroupRef.current);
//         }

//         // -----------------------------------------------------------------------
//         // GENERATE AND ADD PRINCIPAL STRESS AXES
//         // -----------------------------------------------------------------------

//         if (controls.showStressAxes) {
//             const axes = generatePrincipalAxes(
//                 controls.sigmaX,
//                 controls.sigmaY,
//                 controls.sigmaZ
//             );
//             addAxesToScene(axes, sphereGroupRef.current);
//         }
//     }, [
//         controls.sigmaX,
//         controls.sigmaY,
//         controls.sigmaZ,
//         controls.showEquipotential,
//         controls.showSlipLines,
//         controls.showStressAxes,
//         controls.slipLineColor,
//         controls.equipotentialColor
//     ]);

//     // =========================================================================
//     // HELPER FUNCTIONS - THREE.JS GEOMETRY CREATION
//     // =========================================================================

//     /**
//      * Add curve lines to the Three.js scene
//      * Handles both equipotential curves and slip lines
//      */
//     const addLinesToScene = (
//         curves: EquipotentialCurve[] | SlipLineCurve[],
//         color: string,
//         group: THREE.Group
//     ): void => {
//         (curves as any[]).forEach((curve) => {
//             // Extract coordinates from curve
//             const x = (curve as any).x || [];
//             const y = (curve as any).y || [];
//             const z = (curve as any).z || [];

//             if (x.length === 0 || y.length === 0 || z.length === 0) return;
//             if (x.length !== y.length || y.length !== z.length) return;

//             // Create Three.js geometry
//             const geometry = new THREE.BufferGeometry();
//             const positions: number[] = [];

//             for (let i = 0; i < x.length; i++) {
//                 positions.push(x[i], y[i], z[i]);
//             }

//             geometry.setAttribute(
//                 'position',
//                 new THREE.BufferAttribute(new Float32Array(positions), 3)
//             );

//             // Create material and line
//             const material = new THREE.LineBasicMaterial({
//                 color: new THREE.Color(color),
//                 fog: false
//             });

//             const lineObj = new THREE.Line(geometry, material);
//             group.add(lineObj);
//         });
//     };

//     /**
//      * Add principal stress axes to the scene with arrow heads
//      */
//     const addAxesToScene = (axes: PrincipalAxis[], group: THREE.Group): void => {
//         axes.forEach((axis) => {
//             // ---------------------------------------------------------------
//             // MAIN AXIS LINE
//             // ---------------------------------------------------------------

//             const geometry = new THREE.BufferGeometry();
//             const positions = [
//                 0, 0, 0,
//                 axis.axis.x * 1.5, axis.axis.y * 1.5, axis.axis.z * 1.5
//             ];

//             geometry.setAttribute(
//                 'position',
//                 new THREE.BufferAttribute(new Float32Array(positions), 3)
//             );

//             const material = new THREE.LineBasicMaterial({ color: axis.colorHex });
//             const line = new THREE.Line(geometry, material);
//             group.add(line);

//             // ---------------------------------------------------------------
//             // ARROW HEAD
//             // ---------------------------------------------------------------

//             const arrowDir = new THREE.Vector3(axis.axis.x, axis.axis.y, axis.axis.z).normalize();
//             const arrowLength = 0.2;
//             const arrowOrigin = new THREE.Vector3(
//                 axis.axis.x * 1.5,
//                 axis.axis.y * 1.5,
//                 axis.axis.z * 1.5
//             );

//             // Get perpendicular vectors for arrow wings
//             let perpDir1 = new THREE.Vector3()
//                 .crossVectors(arrowDir, new THREE.Vector3(0, 0, 1))
//                 .normalize();

//             if (perpDir1.length() === 0) {
//                 perpDir1.set(1, 0, 0);
//             }

//             const perpDir2 = new THREE.Vector3()
//                 .crossVectors(arrowDir, perpDir1)
//                 .normalize();

//             // Create arrow head geometry
//             const arrowGeometry = new THREE.BufferGeometry();
//             const arrowPositions = [
//                 // First wing
//                 arrowOrigin.x, arrowOrigin.y, arrowOrigin.z,
//                 arrowOrigin.x - arrowDir.x * arrowLength + perpDir1.x * arrowLength * 0.3,
//                 arrowOrigin.y - arrowDir.y * arrowLength + perpDir1.y * arrowLength * 0.3,
//                 arrowOrigin.z - arrowDir.z * arrowLength + perpDir1.z * arrowLength * 0.3,
//                 // Second wing
//                 arrowOrigin.x, arrowOrigin.y, arrowOrigin.z,
//                 arrowOrigin.x - arrowDir.x * arrowLength + perpDir2.x * arrowLength * 0.3,
//                 arrowOrigin.y - arrowDir.y * arrowLength + perpDir2.y * arrowLength * 0.3,
//                 arrowOrigin.z - arrowDir.z * arrowLength + perpDir2.z * arrowLength * 0.3
//             ];

//             arrowGeometry.setAttribute(
//                 'position',
//                 new THREE.BufferAttribute(new Float32Array(arrowPositions), 3)
//             );

//             const arrow = new THREE.LineSegments(arrowGeometry, material);
//             group.add(arrow);
//         });
//     };

//     // =========================================================================
//     // RENDER
//     // =========================================================================

//     return (
//         <div style={{ display: 'flex', height: '100vh', fontFamily: 'Arial, sans-serif' }}>
//             {/* 3D Visualization */}
//             <div
//                 ref={containerRef}
//                 style={{
//                     flex: 1,
//                     position: 'relative'
//                 }}
//             />

//             {/* Control Panel */}
//             <ControlPanel controls={controls} onControlChange={handleControlChange} />
//         </div>
//     );
// };

// // ============================================================================
// // CONTROL PANEL COMPONENT
// // ============================================================================

// interface ControlPanelProps {
//     controls: ControlState;
//     onControlChange: <K extends keyof ControlState>(key: K, value: ControlState[K]) => void;
// }

// const ControlPanel: React.FC<ControlPanelProps> = ({ controls, onControlChange }) => {
//     return (
//         <div
//             style={{
//                 width: 350,
//                 padding: 20,
//                 backgroundColor: '#f5f5f5',
//                 borderLeft: '1px solid #ddd',
//                 overflowY: 'auto'
//             }}
//         >
//             <h2 style={{ marginTop: 0, fontSize: 18, fontWeight: 'bold' }}>
//                 3D Slip Lines & Equipotentials
//             </h2>

//             {/* Stress Tensor Parameters Section */}
//             <div style={{ marginBottom: 25 }}>
//                 <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
//                     Stress Tensor (σ)
//                 </h3>

//                 <StressSlider
//                     label="σ_x (σ₁ - minimum)"
//                     value={controls.sigmaX}
//                     onChange={(val) => onControlChange('sigmaX', val)}
//                 />

//                 <StressSlider
//                     label="σ_y (σ₃ - maximum)"
//                     value={controls.sigmaY}
//                     onChange={(val) => onControlChange('sigmaY', val)}
//                 />

//                 <StressSlider
//                     label="σ_z (σ₂ - intermediate)"
//                     value={controls.sigmaZ}
//                     onChange={(val) => onControlChange('sigmaZ', val)}
//                 />
//             </div>

//             {/* Visualization Options Section */}
//             <div style={{ marginBottom: 25 }}>
//                 <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
//                     Visualization
//                 </h3>

//                 <ToggleControl
//                     label="Show Equipotentials"
//                     checked={controls.showEquipotential}
//                     onChange={(val) => onControlChange('showEquipotential', val)}
//                 >
//                     <ColorPicker
//                         value={controls.equipotentialColor}
//                         onChange={(val) => onControlChange('equipotentialColor', val)}
//                     />
//                 </ToggleControl>

//                 <ToggleControl
//                     label="Show Slip Lines"
//                     checked={controls.showSlipLines}
//                     onChange={(val) => onControlChange('showSlipLines', val)}
//                 >
//                     <ColorPicker
//                         value={controls.slipLineColor}
//                         onChange={(val) => onControlChange('slipLineColor', val)}
//                     />
//                 </ToggleControl>

//                 <ToggleControl
//                     label="Show Stress Axes"
//                     checked={controls.showStressAxes}
//                     onChange={(val) => onControlChange('showStressAxes', val)}
//                 />
//             </div>

//             {/* Info Section */}
//             <div style={{
//                 fontSize: 12,
//                 color: '#666',
//                 borderTop: '1px solid #ddd',
//                 paddingTop: 15
//             }}>
//                 <p style={{ margin: '0 0 8px 0', fontWeight: 'bold' }}>Controls:</p>
//                 <ul style={{ margin: 0, paddingLeft: 16, lineHeight: 1.6 }}>
//                     <li>Drag to rotate</li>
//                     <li>Scroll to zoom</li>
//                 </ul>
//             </div>
//         </div>
//     );
// };

// // ============================================================================
// // UI CONTROL COMPONENTS
// // ============================================================================

// interface StressSliderProps {
//     label: string;
//     value: number;
//     onChange: (value: number) => void;
// }

// const StressSlider: React.FC<StressSliderProps> = ({ label, value, onChange }) => (
//     <div style={{ marginBottom: 15 }}>
//         <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
//             {label}: {value.toFixed(2)}
//         </label>
//         <input
//             type="range"
//             min="-5"
//             max="5"
//             step="0.1"
//             value={value}
//             onChange={(e) => onChange(parseFloat(e.target.value))}
//             style={{ width: '100%', cursor: 'pointer' }}
//         />
//     </div>
// );

// interface ToggleControlProps {
//     label: string;
//     checked: boolean;
//     onChange: (checked: boolean) => void;
//     children?: React.ReactNode;
// }

// const ToggleControl: React.FC<ToggleControlProps> = ({ label, checked, onChange, children }) => (
//     <>
//         <div style={{ marginBottom: 12 }}>
//             <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer', fontSize: 13 }}>
//                 <input
//                     type="checkbox"
//                     checked={checked}
//                     onChange={(e) => onChange(e.target.checked)}
//                     style={{ marginRight: 8, width: 16, height: 16, cursor: 'pointer' }}
//                 />
//                 <span>{label}</span>
//             </label>
//         </div>
//         {checked && children && <div style={{ marginBottom: 12, paddingLeft: 24 }}>{children}</div>}
//     </>
// );

// interface ColorPickerProps {
//     value: string;
//     onChange: (value: string) => void;
// }

// const ColorPicker: React.FC<ColorPickerProps> = ({ value, onChange }) => (
//     <label style={{ display: 'flex', alignItems: 'center', gap: 8, fontSize: 13 }}>
//         Color:
//         <input
//             type="color"
//             value={value}
//             onChange={(e) => onChange(e.target.value)}
//             style={{ width: 40, height: 30, cursor: 'pointer', border: 'none' }}
//         />
//     </label>
// );

// export default SlipLinesVisualization;