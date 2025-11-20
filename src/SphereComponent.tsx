import React, { useEffect, useRef, useState } from 'react';
import * as THREE from 'three';
import { TrackballControls } from 'three/examples/jsm/controls/TrackballControls.js';
import { LineSegments2 } from 'three/examples/jsm/lines/LineSegments2.js';
import { LineGeometry } from 'three/examples/jsm/lines/LineGeometry.js';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';

import {
    generateEquipotentialsBinary,
    generateEquipotentialsParametric,
    generateSpecialEquipotentialMeridian,
    generateSlipLines,
    generatePrincipalAxes,
    getEquipotentialPoint,
    EquipotentialCurve,
    SlipLineCurve,
    PrincipalAxis,
    Point3D,
} from './sphere_utils';
import { BufferGeometry, Float32BufferAttribute, Uint32BufferAttribute } from '../keplerlit/attributes';
import { colorMapNames } from '../keplerlit/colorMap';
import { createIsoContoursFilled, IsoFillReturnedType } from '../keplerlit/IsoContoursFilled';
import { createIsoContourLines } from '../keplerlit/IsoContoursLines';
import { ColorScale } from '../keplerlit/ColorScale';

const computeS2 = (S1: number, S3: number, R: number): number => {
    return S3 + R * (S1 - S3);
};

interface ControlState {
    S1: number;
    S3: number;
    R: number;
    showEquipotential: boolean;
    showSlipLines: boolean;
    showStressAxes: boolean;
    showSpecialEquipotential: boolean;
    slipLineColor: string;
    equipotentialColor: string;
    specialEquipotentialColor: string;
    useParametricEquipotentials: boolean;
    nCurves: number;
    nIsoCurves: number;
    colorTableName: string;
}

const INITIAL_CONTROLS: ControlState = {
    S1: 1.0,
    S3: 3.0,
    R: 0.5,
    showEquipotential: true,
    showSlipLines: true,
    showStressAxes: true,
    showSpecialEquipotential: false,
    slipLineColor: '#D62728',
    equipotentialColor: '#000000',
    specialEquipotentialColor: '#FF00FF',
    useParametricEquipotentials: true,
    nCurves: 12,
    nIsoCurves: 8,
    colorTableName: 'Rainbow'
};

const CAMERA_INITIAL_POS = { x: 3, y: 3, z: 3 };
const SPHERE_OPACITY = 1;
const AXES_SIZE = 0.8;

const SphereComponent: React.FC = () => {
    const containerRef = useRef<HTMLDivElement>(null);
    const sceneRef = useRef<THREE.Scene | null>(null);
    const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
    const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
    const sphereGeometryRef = useRef<THREE.SphereGeometry | null>(null);
    const contourMeshRef = useRef<THREE.Mesh | null>(null);
    const contourLinesRef = useRef<THREE.LineSegments | null>(null);
    const sphereGroupRef = useRef<THREE.Group>(new THREE.Group());
    const colorScaleRef = useRef<ColorScale | null>(null);
    const [controls, setControls] = useState<ControlState>(INITIAL_CONTROLS);

    const handleControlChange = <K extends keyof ControlState>(
        key: K,
        value: ControlState[K]
    ): void => {
        setControls((prev) => ({ ...prev, [key]: value }));
    };

    // =========================================================================
    // SCENE INITIALIZATION
    // =========================================================================

    useEffect(() => {
        if (!containerRef.current) return;

        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0xaaaaaa);
        sceneRef.current = scene;

        const width = containerRef.current.clientWidth;
        const height = containerRef.current.clientHeight;
        const camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        camera.position.set(CAMERA_INITIAL_POS.x, CAMERA_INITIAL_POS.y, CAMERA_INITIAL_POS.z);
        camera.lookAt(0, 0, 0);
        cameraRef.current = camera;

        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(width, height);
        renderer.setPixelRatio(window.devicePixelRatio);
        containerRef.current.appendChild(renderer.domElement);
        rendererRef.current = renderer;

        const ambientLight = new THREE.AmbientLight(0xffffff, 2);
        scene.add(ambientLight);

        const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.99);
        directionalLight1.position.set(5, 5, 5);
        scene.add(directionalLight1);

        const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.99);
        directionalLight2.position.set(-5, -5, -5);
        scene.add(directionalLight2);

        const sphereGeometry = new THREE.SphereGeometry(1, 100, 100);
        sphereGeometryRef.current = sphereGeometry;

        // const background = new GradientBackground(scene, "volcanic");

        // const sphereMaterial = new THREE.MeshPhongMaterial({
        //     color: 0xffffff,
        //     wireframe: false,
        //     opacity: SPHERE_OPACITY,
        //     transparent: true,
        //     vertexColors: true
        // });

        // const sphere = new THREE.Mesh(sphereGeometry, sphereMaterial);
        // scene.add(sphere);

        // const axesHelper = new THREE.AxesHelper(AXES_SIZE);
        // scene.add(axesHelper);

        scene.add(sphereGroupRef.current);

        const trackballControls = new TrackballControls(camera, renderer.domElement);
        trackballControls.rotateSpeed = 1.0;
        trackballControls.zoomSpeed = 1.2;
        trackballControls.panSpeed = 0.8;
        trackballControls.target.set(0, 0, 0);

        // let animationId: number;
        // const animate = (): void => {
        //     animationId = requestAnimationFrame(animate);
        //     trackballControls.update();
        //     renderer.render(scene, camera);
        // };

        let animationId: number;
        let lastFrame = 0;
        const animate = (): void => {
            animationId = requestAnimationFrame(animate);
            const now = Date.now();
            if (now - lastFrame > 33) {  // ~30fps
                trackballControls.update();
                renderer.render(scene, camera);
                lastFrame = now;
            }
        };

        animate();

        const handleResize = (): void => {
            if (!containerRef.current || !cameraRef.current) return;
            const newWidth = containerRef.current.clientWidth;
            const newHeight = containerRef.current.clientHeight;
            cameraRef.current.aspect = newWidth / newHeight;
            cameraRef.current.updateProjectionMatrix();
            renderer.setSize(newWidth, newHeight);

            // Update color scale position and size on resize
            if (colorScaleRef.current) {
                const x = 50;
                const y = 50;
                colorScaleRef.current.updatePosition(x, y);
                colorScaleRef.current.updateSize(30, 200);
            }
        };

        window.addEventListener('resize', handleResize);

        // Initialize color scale
        const colorStops = [
            { position: 0, color: '#0066cc' },
            { position: 0.5, color: '#ffff00' },
            { position: 1, color: '#cc0000' }
        ]

        const colorScale = new ColorScale({
            canvas: renderer.domElement,
            x: 50, //window.innerWidth - 100,
            y: 50,
            width: 30,
            height: 200,
            min: 0,
            max: 100,
            attributeName: 'Stress',
            orientation: 'vertical',
            colorMapName: INITIAL_CONTROLS.colorTableName,
            autoRender: true,
            colorStops
        });
        colorScaleRef.current = colorScale;

        return () => {
            window.removeEventListener('resize', handleResize);
            trackballControls.dispose();
            cancelAnimationFrame(animationId);
            renderer.dispose();
            if (containerRef.current?.contains(renderer.domElement)) {
                containerRef.current.removeChild(renderer.domElement);
            }
        };
    }, []);

    // =========================================================================
    // VISUALIZATION UPDATE
    // =========================================================================

    const curveResolution = 200;
    useEffect(() => {
        if (!sphereGroupRef.current) return;
        sphereGroupRef.current.clear();

        if (controls.showEquipotential) {
            const config = {
                nCurves: controls.nCurves,
                nIsoCurves: controls.nIsoCurves,
                curveResolution
            };
            const S2 = computeS2(controls.S1, controls.S3, controls.R);
            let equipotentials = controls.useParametricEquipotentials
                ? generateEquipotentialsParametric(controls.S1, controls.S3, S2, config)
                : generateEquipotentialsBinary(controls.S1, controls.S3, S2, config);
            addLinesToScene(equipotentials, controls.equipotentialColor, sphereGroupRef.current, false);
        }

        if (controls.showSpecialEquipotential) {
            const S2 = computeS2(controls.S1, controls.S3, controls.R);
            const config = { curveResolution: 500 };
            const specialMeridians = generateSpecialEquipotentialMeridian(controls.S1, controls.S3, S2, config);
            addLinesToScene(specialMeridians, controls.specialEquipotentialColor, sphereGroupRef.current, false, true);
        }

        if (controls.showSlipLines) {
            const config = {
                nCurves: controls.nCurves,
                nIsoCurves: controls.nIsoCurves,
                curveResolution
            };
            const S2 = computeS2(controls.S1, controls.S3, controls.R);
            const slipLines = generateSlipLines(
                controls.S1, controls.S3, S2, config
            );
            addLinesToScene(slipLines, controls.slipLineColor, sphereGroupRef.current, true);
        }

        if (controls.showStressAxes) {
            const S2 = computeS2(controls.S1, controls.S3, controls.R);
            const axes = generatePrincipalAxes(controls.S1, controls.S3, S2);
            addAxesToScene(axes, sphereGroupRef.current);
        }

    }, [
        controls.S1,
        controls.S3,
        controls.R,
        controls.showEquipotential,
        controls.showSlipLines,
        controls.showStressAxes,
        controls.slipLineColor,
        controls.equipotentialColor,
        controls.nCurves,
        controls.nIsoCurves,
        controls.colorTableName
        // colorScaleRef.current
    ]);

    // EFFECT 1B: Update color scale when color table changes
    useEffect(() => {
        if (colorScaleRef.current && controls.colorTableName) {
            colorScaleRef.current.setColorMap(controls.colorTableName);
        }
    }, [controls.colorTableName, colorScaleRef.current]);

    // EFFECT 2: Update sphere coloring ONLY (depends on stress values only)
    useEffect(() => {
        if (!sphereGeometryRef.current) return;

        const geom = sphereGeometryRef.current;
        const position = geom.attributes.position;

        /*
        let attribute: THREE.BufferAttribute | null = geom.getAttribute('equipotentialValue') as THREE.BufferAttribute;
        if (!attribute) {
            const equipotentialValues = new Float32Array(position.count);
            geom.setAttribute('equipotentialValue', new THREE.BufferAttribute(equipotentialValues, 1));
            attribute = geom.getAttribute('equipotentialValue') as THREE.BufferAttribute;
        }

        const attr: number[] = [];
        const S2 = computeS2(controls.S1, controls.S3, controls.R);

        for (let i = 0; i < position.count; i++) {
            const x = position.getX(i);
            const y = position.getY(i);
            const z = position.getZ(i);
            const v = getEquipotentialPoint({ x, y, z }, controls.S1, S2, controls.S3);
            attr.push(v);
        }

        const mm = minMax(attr);
        const normalizeAttr = (v: number) => (v - mm[0]) / (mm[1] - mm[0]);
        const lutTable = createLut("Igeoss", 512);
        lutTable.setMin(mm[0]);
        lutTable.setMax(mm[1]);

        const colors: number[] = [];
        attr.forEach((v) => {
            const c = fromValueToColor(normalizeAttr(v), {
                min: mm[0],
                max: mm[1],
                defaultColor: new Color("#ffffff"),
                lutTable
            });
            colors.push(...c);
        });

        geom.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        geom.getAttribute("color").needsUpdate = true;
        */

        const generateIndices = (vertexCount: number) => {
            const indices = [];
            for (let i = 0; i < vertexCount - 2; i += 3) {
                indices.push(i, i + 1, i + 2);
            }
            return indices;
        }

        const convertToKeplerGeometry = (threeGeometry: THREE.BufferGeometry) => {
            const positions = threeGeometry.attributes.position.array;
            const indices = threeGeometry.index ? threeGeometry.index.array : generateIndices(positions.length / 3);

            const keplerPositions = new Float32BufferAttribute(Array.from(positions), 3);
            const keplerIndices = new Uint32BufferAttribute(Array.from(indices), 1);

            const keplerGeometry = new BufferGeometry();
            keplerGeometry.setPositions(keplerPositions);
            keplerGeometry.setIndices(keplerIndices);

            return keplerGeometry;
        }

        const createMeshFromKeplerResult = (result: IsoFillReturnedType) => {
            const geometry = new THREE.BufferGeometry();

            // Convert positions
            geometry.setAttribute('position', new THREE.Float32BufferAttribute(result.position, 3));

            // Convert indices
            geometry.setIndex(new THREE.Uint32BufferAttribute(result.index, 1));

            // Convert colors
            const colors = new Float32Array(result.color);
            geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

            geometry.computeVertexNormals();
            // geometry.setAttribute('normal', new THREE.Float32BufferAttribute(result.normal, 3))

            const material = new THREE.MeshPhongMaterial({
                vertexColors: true,
                side: THREE.DoubleSide,
                wireframe: false,
                flatShading: false,
                polygonOffset: true,
                polygonOffsetFactor: .5
            });

            return new THREE.Mesh(geometry, material);
        }

        // Generate scalar field based on height and a wave pattern
        const scalarField = [];
        for (let i = 0; i < position.count; ++i) {
            const point = {
                x: position.getX(i),
                y: position.getY(i),
                z: position.getZ(i)
            }

            const S1 = controls.S1;
            const S3 = controls.S3;
            const S2 = computeS2(S1, S3, controls.R);

            // Scalar field combining height and radial distance
            const scalar = getEquipotentialPoint(point, S3, S2, S1)
            scalarField.push(scalar);
        }

        const numContours = controls.nIsoCurves + 1;

        const minVal = Math.min(...scalarField);
        const maxVal = Math.max(...scalarField);
        const isoList = [];
        for (let i = 0; i < numContours; i++) {
            isoList.push(minVal + (i / (numContours - 1)) * (maxVal - minVal));
        }

        colorScaleRef.current?.updateRange(minVal, maxVal);

        // Convert Three.js geometry to kepler format
        const keplerGeometry = convertToKeplerGeometry(geom);

        const isos = createIsoContoursFilled(keplerGeometry, scalarField, isoList, {
            min: minVal,
            max: maxVal,
            lut: colorScaleRef.current?.getColorMapName() || "Rainbow",
            nbColors: 128
        })

        const mesh = createMeshFromKeplerResult(isos!);
        if (contourMeshRef.current) {
            sphereGroupRef.current.remove(contourMeshRef.current);
        }
        contourMeshRef.current = mesh;
        sphereGroupRef.current.add(mesh);

        // Iso lines
        if (0) {
            const isoLines = createIsoContourLines(keplerGeometry, scalarField, isoList, "#000000", "Rainbow");
            if (isoLines.positions.length > 0) {
                const lineGeometry = new THREE.BufferGeometry();
                lineGeometry.setAttribute('position', new THREE.Float32BufferAttribute(isoLines.positions, 3));

                // const colors = new Float32Array(result.color);
                // lineGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

                const lineMaterial = new THREE.LineBasicMaterial({
                    // vertexColors: true,
                    color: 0x000000,
                    linewidth: 2
                });

                const contourLines = new THREE.LineSegments(lineGeometry, lineMaterial);
                // this.scene.add(this.contourLines);
                if (contourLinesRef.current) {
                    sphereGroupRef.current.remove(contourLinesRef.current);
                }
                contourLinesRef.current = contourLines;
                sphereGroupRef.current.add(contourLines);
            }
        }

    }, [
        controls.S1,
        controls.S3,
        controls.R,
        controls.showEquipotential,
        controls.showSlipLines,
        controls.showStressAxes,
        controls.slipLineColor,
        controls.equipotentialColor,
        controls.nCurves,
        controls.nIsoCurves,
        controls.colorTableName
        // colorScaleRef.current
    ]);  // â† ONLY stress values!

    // =========================================================================
    // HELPER FUNCTIONS
    // =========================================================================

    const addLinesToScene = (
        curves: EquipotentialCurve[] | SlipLineCurve[],
        color: string,
        group: THREE.Group,
        isSlipLine: boolean,
        isSpecial: boolean = false,
        lineWidth = 2
    ): void => {
        let renderedCount = 0;

        (curves as any[]).forEach((curve) => {
            const x = (curve as any).x || [];
            const y = (curve as any).y || [];
            const z = (curve as any).z || [];

            if (x.length === 0 || y.length === 0 || z.length === 0) return;
            if (x.length !== y.length || y.length !== z.length) {
                console.warn(`Curve: Length mismatch`);
                return;
            }

            const positions: number[] = [];
            const scale = 1.001;  // Offset to avoid Z-fighting with sphere

            for (let i = 0; i < x.length; i++) {
                if (!isFinite(x[i]) || !isFinite(y[i]) || !isFinite(z[i])) {
                    return;
                }

                // Three.js expects: x = horizontal, y = vertical, z = horizontal
                // Swap (x,y,z) to (y,z,x) equivalent to (North, Up, East) for display
                positions.push(y[i] * scale, z[i] * scale, x[i] * scale);
            }

            if (positions.length === 0) return;

            const linewidth = lineWidth//isSpecial ? 10 : (isSlipLine ? 6 : 4);

            const lineGeometry = new LineGeometry().setPositions(positions);
            const lineMaterial = new LineMaterial({
                color: new THREE.Color(color),
                linewidth: linewidth,
                resolution: new THREE.Vector2(
                    rendererRef.current!.domElement.clientWidth,
                    rendererRef.current!.domElement.clientHeight
                )
            });
            const lineObj = new LineSegments2(lineGeometry, lineMaterial);
            group.add(lineObj);

            renderedCount++;
        });

        // console.log(`Rendered ${renderedCount} curves`);
    };

    const addAxesToScene = (axes: PrincipalAxis[], group: THREE.Group): void => {
        axes.forEach((axis) => {
            const geometry = new THREE.BufferGeometry();
            const positions = [
                0, 0, 0,
                axis.axis.y * 1.5,   // North
                axis.axis.z * 1.5,   // Up
                axis.axis.x * 1.5    // East
            ];

            geometry.setAttribute(
                'position',
                new THREE.BufferAttribute(new Float32Array(positions), 3)
            );

            const material = new THREE.LineBasicMaterial({
                color: axis.colorHex,
                polygonOffset: true,
                polygonOffsetFactor: 1,
                polygonOffsetUnits: 1
            });
            const line = new THREE.Line(geometry, material);
            group.add(line);

            const arrowDir = new THREE.Vector3(axis.axis.y, axis.axis.z, axis.axis.x).normalize();
            const arrowLength = 0.2;
            const arrowOrigin = new THREE.Vector3(
                axis.axis.y * 1.5,
                axis.axis.z * 1.5,
                axis.axis.x * 1.5
            );

            let perpDir1 = new THREE.Vector3()
                .crossVectors(arrowDir, new THREE.Vector3(0, 0, 1))
                .normalize();

            if (perpDir1.length() === 0) {
                perpDir1.set(1, 0, 0);
            }

            const perpDir2 = new THREE.Vector3()
                .crossVectors(arrowDir, perpDir1)
                .normalize();

            const arrowGeometry = new THREE.BufferGeometry();
            const arrowPositions = [
                arrowOrigin.x, arrowOrigin.y, arrowOrigin.z,
                arrowOrigin.x - arrowDir.x * arrowLength + perpDir1.x * arrowLength * 0.3,
                arrowOrigin.y - arrowDir.y * arrowLength + perpDir1.y * arrowLength * 0.3,
                arrowOrigin.z - arrowDir.z * arrowLength + perpDir1.z * arrowLength * 0.3,
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

    return (
        <div style={{ display: 'flex', height: '100vh', fontFamily: 'Arial, sans-serif' }}>
            <div
                ref={containerRef}
                style={{
                    flex: 1,
                    position: 'relative'
                }}
            />
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
                width: 380,
                padding: 20,
                backgroundColor: '#f5f5f5',
                borderLeft: '1px solid #ddd',
                overflowY: 'auto'
            }}
        >
            <h2 style={{ marginTop: 0, fontSize: 18, fontWeight: 'bold' }}>
                3D Slip Lines & Equipotentials
            </h2>

            {/* ================================================================ */}
            {/* PRINCIPAL STRESSES SECTION */}
            {/* ================================================================ */}

            <div style={{ marginBottom: 25 }}>
                <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
                    Principal Stresses
                </h3>

                <StressSlider
                    label="Sâ‚ (Ïƒâ‚ - minimum)"
                    value={controls.S1}
                    onChange={(val) => onControlChange('S1', val)}
                    min={-5}
                    max={5}
                    step={0.1}
                />

                <StressSlider
                    label="Sâ‚ƒ (Ïƒâ‚ƒ - maximum)"
                    value={controls.S3}
                    onChange={(val) => onControlChange('S3', val)}
                    min={-5}
                    max={5}
                    step={0.1}
                />

                <div style={{ marginBottom: 15 }}>
                    <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
                        R (stress ratio): {controls.R.toFixed(2)}
                    </label>
                    <input
                        type="range"
                        min="0"
                        max="1"
                        step="0.01"
                        value={controls.R}
                        onChange={(e) => onControlChange('R', parseFloat(e.target.value))}
                        style={{ width: '100%', cursor: 'pointer' }}
                    />
                    <small style={{ color: '#999', fontSize: 11 }}>
                        R = (Sâ‚‚ - Sâ‚ƒ) / (Sâ‚ - Sâ‚ƒ), where Sâ‚‚ = {computeS2(controls.S1, controls.S3, controls.R).toFixed(3)}
                    </small>
                </div>
            </div>

            {/* ================================================================ */}
            {/* COMPUTATION PARAMETERS SECTION */}
            {/* ================================================================ */}

            {/* <div style={{ marginBottom: 25 }}>
                <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
                    Computation Settings
                </h3>

                <div style={{ marginBottom: 15 }}>
                    <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer', fontSize: 13, marginBottom: 8 }}>
                        <input
                            type="checkbox"
                            checked={controls.useParametricEquipotentials}
                            onChange={(e) => onControlChange('useParametricEquipotentials', e.target.checked)}
                            style={{ marginRight: 8, width: 16, height: 16, cursor: 'pointer' }}
                        />
                        <span>Use Parametric Equipotentials</span>
                    </label>
                    <small style={{ color: '#999', fontSize: 11, marginLeft: 24 }}>
                        {controls.useParametricEquipotentials
                            ? 'Improved method using Mohr circle parametrization'
                            : 'Binary search method (stable baseline)'}
                    </small>
                </div>

                <div style={{ marginBottom: 15 }}>
                    <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
                        # Slip Lines: {controls.nCurves}
                    </label>
                    <input
                        type="range"
                        min="4"
                        max="20"
                        step="1"
                        value={controls.nCurves}
                        onChange={(e) => onControlChange('nCurves', parseInt(e.target.value))}
                        style={{ width: '100%', cursor: 'pointer' }}
                    />
                    <small style={{ color: '#999', fontSize: 11 }}>More curves = denser visualization</small>
                </div>

                <div style={{ marginBottom: 15 }}>
                    <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
                        # Equipotentials: {controls.nIsoCurves}
                    </label>
                    <input
                        type="range"
                        min="4"
                        max="30"
                        step="1"
                        value={controls.nIsoCurves}
                        onChange={(e) => onControlChange('nIsoCurves', parseInt(e.target.value))}
                        style={{ width: '100%', cursor: 'pointer' }}
                    />
                    <small style={{ color: '#999', fontSize: 11 }}>More levels = finer stress detail</small>
                </div>
            </div> */}

            {/* ================================================================ */}
            {/* VISUALIZATION TOGGLES SECTION */}
            {/* ================================================================ */}

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
                    <div style={{ marginBottom: 15 }}>
                        <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
                            # Equipotentials: {controls.nIsoCurves}
                        </label>
                        <input
                            type="range"
                            min="4"
                            max="30"
                            step="1"
                            value={controls.nIsoCurves}
                            onChange={(e) => onControlChange('nIsoCurves', parseInt(e.target.value))}
                            style={{ width: '100%', cursor: 'pointer' }}
                        />
                        <small style={{ color: '#999', fontSize: 11 }}>More levels = finer stress detail</small>
                    </div>
                </ToggleControl>

                <ToggleControl
                    label="Show Special Equipotential"
                    checked={controls.showSpecialEquipotential}
                    onChange={(val) => onControlChange('showSpecialEquipotential', val)}
                >
                    <ColorPicker
                        value={controls.specialEquipotentialColor}
                        onChange={(val) => onControlChange('specialEquipotentialColor', val)}
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
                    <div style={{ marginBottom: 15 }}>
                        <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
                            # Slip Lines: {controls.nCurves}
                        </label>
                        <input
                            type="range"
                            min="4"
                            max="20"
                            step="1"
                            value={controls.nCurves}
                            onChange={(e) => onControlChange('nCurves', parseInt(e.target.value))}
                            style={{ width: '100%', cursor: 'pointer' }}
                        />
                        <small style={{ color: '#999', fontSize: 11 }}>More curves = denser visualization</small>
                    </div>
                </ToggleControl>

                <ToggleControl
                    label="Show Stress Axes"
                    checked={controls.showStressAxes}
                    onChange={(val) => onControlChange('showStressAxes', val)}
                />
            </div>

            {/* ================================================================ */}
            {/* COLOR SCALE SECTION */}
            {/* ================================================================ */}

            <div style={{ marginBottom: 25 }}>
                <h3 style={{ fontSize: 14, fontWeight: 'bold', marginBottom: 12 }}>
                    Color Scale
                </h3>

                <ComboControl
                    label="Color Table"
                    value={controls.colorTableName}
                    options={colorMapNames()}
                    onChange={(val) => onControlChange('colorTableName', val)}
                />
            </div>

            {/* ================================================================ */}
            {/* CONTROLS INFO SECTION */}
            {/* ================================================================ */}

            <div style={{
                fontSize: 12,
                color: '#666',
                borderTop: '1px solid #ddd',
                paddingTop: 15
            }}>
                <p style={{ margin: '0 0 8px 0', fontWeight: 'bold' }}>Mouse Controls:</p>
                <ul style={{ margin: 0, paddingLeft: 16, lineHeight: 1.6 }}>
                    <li>Left drag: Rotate</li>
                    <li>Middle drag: Zoom</li>
                    <li>Right drag: Pan</li>
                    <li>Scroll: Zoom</li>
                </ul>
                <p style={{ margin: '12px 0 8px 0', fontSize: 11, color: '#999', fontStyle: 'italic' }}>
                    Open browser console (F12) to see computation logs
                </p>
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
    min?: string | number;
    max?: string | number;
    step?: string | number;
}

const StressSlider: React.FC<StressSliderProps> = ({
    label,
    value,
    onChange,
    min = "-5",
    max = "5",
    step = "0.1"
}) => (
    <div style={{ marginBottom: 15 }}>
        <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
            {label}: {value.toFixed(2)}
        </label>
        <input
            type="range"
            min={min}
            max={max}
            step={step}
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

interface ComboControlProps {
    label: string;
    value: string;
    options: string[];
    onChange: (value: string) => void;
}

const ComboControl: React.FC<ComboControlProps> = ({ label, value, options, onChange }) => (
    <div style={{ marginBottom: 15 }}>
        <label style={{ display: 'block', marginBottom: 5, fontSize: 13 }}>
            {label}
        </label>
        <select
            value={value}
            onChange={(e) => onChange(e.target.value)}
            style={{
                width: '100%',
                padding: '6px 8px',
                fontSize: 13,
                cursor: 'pointer',
                border: '1px solid #ccc',
                borderRadius: '3px',
                backgroundColor: '#fff'
            }}
        >
            {options.map((option) => (
                <option key={option} value={option}>
                    {option}
                </option>
            ))}
        </select>
    </div>
);

export default SphereComponent;