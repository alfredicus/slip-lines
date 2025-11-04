import * as THREE from 'three';
import { Text } from 'troika-three-text';

/* USAGE
    import * as THREE from 'three';
    import { BoundingBoxAxesTroika } from './BoundingBoxAxesTroika';

    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(60, innerWidth/innerHeight, 0.01, 100);
    camera.position.set(2.5, 1.8, 2.6);

    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(innerWidth, innerHeight);
    document.body.appendChild(renderer.domElement);

    // Demo model
    const model = new THREE.Mesh(
        new THREE.TorusKnotGeometry(0.5, 0.18, 200, 24),
        new THREE.MeshStandardMaterial({ color: 0x66ccff })
    );
    scene.add(model);

    // Helper
    const helper = new BoundingBoxAxesTroika(model, {
        origin: 'min',              // start at min corner
        billboardLabels: true,      // face labels toward the camera
        addAxisEndLabels: true,
        tickTarget: 10,
        fontSize: 0.06
    });
    scene.add(helper);

    const light = new THREE.DirectionalLight(0xffffff, 1.0);
    light.position.set(3, 4, 2);
    scene.add(light);

    function onResize() {
        camera.aspect = innerWidth / innerHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(innerWidth, innerHeight);
    }
    addEventListener('resize', onResize);

    renderer.setAnimationLoop(() => {
        model.rotation.y += 0.005;
        helper.update(camera); // keeps box in sync + billboarding
        renderer.render(scene, camera);
    });
*/

export type AxisName = 'X' | 'Y' | 'Z';

export interface BoundingBoxAxesOptions {
    tickTarget?: number;                 // ~desired ticks per edge (we “nice” it)
    tickSizeRatio?: number;              // tick length relative to longest box side
    fontSize?: number;                   // world units
    colors?: { x: number; y: number; z: number; box: number; text: number };
    edges?: AxisName[];                  // which edges to draw (default: X,Y,Z)
    origin?: 'min' | 'max';              // corner to start from
    showBox?: boolean;                   // show wireframe box
    labelOffset?: number;                // offset from edge (world units)
    addAxisEndLabels?: boolean;          // add "X"/"Y"/"Z" at edge ends
    billboardLabels?: boolean;           // make labels face camera on update()
    format?: (value: number, step: number) => string; // custom tick formatter
}

export class BoundingBoxAxesTroika extends THREE.Group {
    private target: THREE.Object3D;
    private readonly opts: Required<BoundingBoxAxesOptions>;
    private boxHelper?: THREE.Box3Helper;
    private box: THREE.Box3 = new THREE.Box3();
    private size: THREE.Vector3 = new THREE.Vector3();
    private maxSide = 1;

    constructor(target: THREE.Object3D, options: BoundingBoxAxesOptions = {}) {
        super();
        this.target = target;

        this.opts = {
            tickTarget: options.tickTarget ?? 10,
            tickSizeRatio: options.tickSizeRatio ?? 0.02,
            fontSize: options.fontSize ?? 0.05,
            colors: options.colors ?? {
                x: 0xff5555,
                y: 0x55ff55,
                z: 0x5599ff,
                box: 0x888888,
                text: 0xffffff
            },
            edges: options.edges ?? ['X', 'Y', 'Z'],
            origin: options.origin ?? 'min',
            showBox: options.showBox ?? true,
            labelOffset: options.labelOffset ?? 0.06,
            addAxisEndLabels: options.addAxisEndLabels ?? true,
            billboardLabels: options.billboardLabels ?? false,
            format: options.format ?? ((v, step) => {
                // Reasonable default formatter
                const dec = Math.max(0, -Math.floor(Math.log10(step)));
                return v.toFixed(Math.min(dec, 4));
            })
        };

        this.rebuild();
    }

    setTarget(target: THREE.Object3D) {
        this.target = target;
        this.rebuild();
    }

    /** Rebuild the helper from the current target’s bounds. Call if geometry changes drastically. */
    public rebuild(): void {
        this.clearNodes();

        this.box.setFromObject(this.target);
        if (!isFinite(this.box.min.x) || !isFinite(this.box.max.x)) return;

        this.box.getSize(this.size);
        this.maxSide = Math.max(this.size.x, this.size.y, this.size.z);

        if (this.opts.showBox) {
            this.boxHelper = new THREE.Box3Helper(this.box.clone(), this.opts.colors.box);
            this.add(this.boxHelper);
        }

        const originCorner = (this.opts.origin === 'min') ? this.box.min.clone() : this.box.max.clone();

        // For edges we want three, mutually orthogonal directions from corner
        const dirs: Record<AxisName, THREE.Vector3> = {
            X: (this.opts.origin === 'min') ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(-1, 0, 0),
            Y: (this.opts.origin === 'min') ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(0, -1, 0),
            Z: (this.opts.origin === 'min') ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(0, 0, -1)
        };
        const lengths: Record<AxisName, number> = { X: this.size.x, Y: this.size.y, Z: this.size.z };
        const axisColors: Record<AxisName, number> = {
            X: this.opts.colors.x, Y: this.opts.colors.y, Z: this.opts.colors.z
        };

        // Outward directions to push ticks/labels "outside" the box (away from the cube)
        // If origin=min we push toward -orthos; if origin=max we push toward +orthos.
        const outwardSign = (this.opts.origin === 'min') ? -1 : +1;

        // Build requested edges
        for (const axis of this.opts.edges) {
            const edgeGroup = new THREE.Group();
            this.add(edgeGroup);

            const dir = dirs[axis].clone();
            const length = lengths[axis];
            const color = axisColors[axis];

            // Edge line
            const a = originCorner.clone();
            const b = originCorner.clone().add(dir.clone().multiplyScalar(length));
            edgeGroup.add(this.makeLine(a, b, color));

            // Orthogonal vectors for tick strokes (choose the two axes different from current)
            const ortho1 = (axis === 'X') ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
            const ortho2 = (axis === 'Z') ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 0, 1);

            const tickLen = this.opts.tickSizeRatio * this.maxSide;
            const step = this.niceStep(length / this.opts.tickTarget);
            const start = 0;
            const end = length;
            const firstTick = Math.ceil(start / step) * step;

            for (let t = firstTick; t <= end + 1e-6; t += step) {
                const p = a.clone().add(dir.clone().multiplyScalar(t));

                // two short strokes
                const s1a = p.clone().add(ortho1.clone().multiplyScalar(+tickLen * 0.5 * outwardSign));
                const s1b = p.clone().add(ortho1.clone().multiplyScalar(-tickLen * 0.5 * outwardSign));
                const s2a = p.clone().add(ortho2.clone().multiplyScalar(+tickLen * 0.5 * outwardSign));
                const s2b = p.clone().add(ortho2.clone().multiplyScalar(-tickLen * 0.5 * outwardSign));
                edgeGroup.add(this.makeLine(s1a, s1b, color));
                edgeGroup.add(this.makeLine(s2a, s2b, color));

                // numeric label (slightly outside along both orthos)
                const labelPos = p.clone()
                    .add(ortho1.clone().multiplyScalar(this.opts.labelOffset * outwardSign))
                    .add(ortho2.clone().multiplyScalar(this.opts.labelOffset * 0.6 * outwardSign));
                const text = this.opts.format(t, step);
                edgeGroup.add(this.makeText(text, labelPos, this.opts.fontSize, this.opts.colors.text));
            }

            // Axis end label (e.g., "X")
            if (this.opts.addAxisEndLabels) {
                const endPos = b.clone()
                    .add(ortho1.clone().multiplyScalar(this.opts.labelOffset * 0.6 * outwardSign))
                    .add(ortho2.clone().multiplyScalar(this.opts.labelOffset * 0.6 * outwardSign));
                edgeGroup.add(this.makeText(axis, endPos, this.opts.fontSize * 1.3, color));
            }
        }
    }

    /** Call every frame if billboardLabels=true, and/or when the target transforms change. */
    public update(camera?: THREE.Camera): void {
        // Keep the wire box in sync
        this.box.setFromObject(this.target);
        if (this.boxHelper) this.boxHelper.box.copy(this.box);

        // Billboard labels if requested
        if (this.opts.billboardLabels && camera) {
            const camPos = new THREE.Vector3().setFromMatrixPosition(camera.matrixWorld);
            this.traverse(obj => {
                if ((obj as any).isTroikaText) {
                    const m = obj as Text;
                    // Simple billboard: face camera position
                    obj.lookAt(camPos);
                    // Troika needs sync if text or params changed; rotation is fine without sync.
                }
            });
        }
    }

    /** Clean up GPU resources. */
    public dispose(): void {
        this.traverse(obj => {
            if (obj instanceof THREE.Line) {
                obj.geometry.dispose();
                (obj.material as THREE.Material)?.dispose?.();
            }
            if ((obj as any).isTroikaText) {
                (obj as Text).dispose();
            }
        });
        this.clear();
    }

    // ---------- helpers ----------

    private makeLine(a: THREE.Vector3, b: THREE.Vector3, color: number): THREE.Line {
        const geom = new THREE.BufferGeometry().setFromPoints([a, b]);
        const mat = new THREE.LineBasicMaterial({ color });
        return new THREE.Line(geom, mat);
    }

    private makeText(text: string, pos: THREE.Vector3, fontSize: number, color: number): Text {
        const t = new Text();
        (t as any).isTroikaText = true; // tag for update()
        t.text = text;
        t.fontSize = fontSize;
        t.color = color;
        t.anchorX = 'center';
        t.anchorY = 'middle';
        t.position.copy(pos);
        t.depthTest = true;   // occlude behind geometry
        t.depthWrite = true;
        t.sync();             // build geometry
        return t;
    }

    private niceStep(step: number): number {
        const pow10 = Math.pow(10, Math.floor(Math.log10(step)));
        const n = step / pow10;
        const nice = (n < 1.5) ? 1 : (n < 3.5) ? 2 : (n < 7.5) ? 5 : 10;
        return nice * pow10;
    }

    private clearNodes(): void {
        // Dispose of any previous geometry/materials/texts to avoid leaks
        this.traverse(obj => {
            if (obj instanceof THREE.Line) {
                obj.geometry.dispose();
                (obj.material as THREE.Material)?.dispose?.();
            }
            if ((obj as any).isTroikaText) {
                (obj as Text).dispose();
            }
        });
        // Remove all children from the group
        this.clear();
    }
}
