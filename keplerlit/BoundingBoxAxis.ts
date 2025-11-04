import * as THREE from 'three';
import { TextGeometry } from 'three/examples/jsm/geometries/TextGeometry.js';
import { FontLoader } from 'three/examples/jsm/loaders/FontLoader.js';

/* USAGE
    import * as THREE from 'three';
    import { BoundingBoxAxes } from './BoundingBoxAxes';

    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(60, innerWidth / innerHeight, 0.01, 100);
    camera.position.set(2, 1.5, 2);
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(innerWidth, innerHeight);
    document.body.appendChild(renderer.domElement);

    const model = new THREE.Mesh(
        new THREE.TorusKnotGeometry(0.4, 0.12, 120, 16),
        new THREE.MeshStandardMaterial({ color: 0x66ccff })
    );
    scene.add(model);

    const helper = new BoundingBoxAxes(model);
    scene.add(helper);

    const light = new THREE.DirectionalLight(0xffffff, 1);
    light.position.set(2, 2, 1);
    scene.add(light);

    function animate() {
        requestAnimationFrame(animate);
        model.rotation.y += 0.01;
        helper.update();
        renderer.render(scene, camera);
    }
    animate();
*/

/**
 * Create a bounding box around an object, with labeled, graduated axes.
 * Pure Three.js version — no CSS2DRenderer or external text libs.
 */
export class BoundingBoxAxes extends THREE.Group {
    private boxHelper: THREE.Box3Helper;
    private font: THREE.Font | null = null;
    private params: {
        tickCount: number;
        fontSize: number;
        tickLength: number;
        colors: { x: number; y: number; z: number; box: number };
    };

    constructor(
        private target: THREE.Object3D,
        fontUrl = 'https://threejs.org/examples/fonts/helvetiker_regular.typeface.json',
        params?: Partial<{
            tickCount: number;
            fontSize: number;
            tickLength: number;
            colors: { x: number; y: number; z: number; box: number };
        }>
    ) {
        super();

        this.params = {
            tickCount: params?.tickCount ?? 10,
            fontSize: params?.fontSize ?? 0.05,
            tickLength: params?.tickLength ?? 0.02,
            colors: params?.colors ?? {
                x: 0xff5555,
                y: 0x55ff55,
                z: 0x5599ff,
                box: 0x888888,
            },
        };

        this.boxHelper = new THREE.Box3Helper(
            new THREE.Box3().setFromObject(target),
            this.params.colors.box
        );
        this.add(this.boxHelper);

        new FontLoader().load(fontUrl, (font) => {
            this.font = font;
            this.build();
        });
    }

    private build() {
        const box = new THREE.Box3().setFromObject(this.target);
        const size = new THREE.Vector3();
        const min = box.min.clone();
        const max = box.max.clone();
        box.getSize(size);

        // Edges origin corner (min corner)
        this.addEdge(min, new THREE.Vector3(1, 0, 0), size.x, 'X', this.params.colors.x);
        this.addEdge(min, new THREE.Vector3(0, 1, 0), size.y, 'Y', this.params.colors.y);
        this.addEdge(min, new THREE.Vector3(0, 0, 1), size.z, 'Z', this.params.colors.z);
    }

    private addEdge(
        origin: THREE.Vector3,
        dir: THREE.Vector3,
        length: number,
        axisName: 'X' | 'Y' | 'Z',
        color: number
    ) {
        const group = new THREE.Group();
        this.add(group);

        // Edge line
        const geom = new THREE.BufferGeometry().setFromPoints([
            origin,
            origin.clone().add(dir.clone().multiplyScalar(length)),
        ]);
        const line = new THREE.Line(geom, new THREE.LineBasicMaterial({ color }));
        group.add(line);

        const tickCount = this.params.tickCount;
        const tickLen = this.params.tickLength;
        const tickStep = this.niceStep(length / tickCount);

        // Orthogonal directions for small ticks
        const ortho1 =
            axisName === 'X'
                ? new THREE.Vector3(0, 1, 0)
                : axisName === 'Y'
                    ? new THREE.Vector3(1, 0, 0)
                    : new THREE.Vector3(0, 1, 0);
        const ortho2 =
            axisName === 'Z'
                ? new THREE.Vector3(1, 0, 0)
                : new THREE.Vector3(0, 0, 1);

        for (let t = 0; t <= length + 1e-6; t += tickStep) {
            const p = origin.clone().add(dir.clone().multiplyScalar(t));
            const a = p.clone().add(ortho1.clone().multiplyScalar(tickLen / 2));
            const b = p.clone().add(ortho1.clone().multiplyScalar(-tickLen / 2));
            const tickGeom = new THREE.BufferGeometry().setFromPoints([a, b]);
            group.add(new THREE.Line(tickGeom, new THREE.LineBasicMaterial({ color })));

            this.addLabel(group, p, t, axisName, color);
        }

        // Axis end label
        const endPos = origin.clone().add(dir.clone().multiplyScalar(length + tickLen * 2));
        this.addText(group, axisName, endPos, color);
    }

    private addLabel(
        group: THREE.Group,
        pos: THREE.Vector3,
        value: number,
        axis: string,
        color: number
    ) {
        const str = this.formatNumber(value);
        const offset = 0.04;
        const labelPos = pos.clone().add(new THREE.Vector3(0, offset, 0));
        this.addText(group, str, labelPos, color);
    }

    private addText(
        group: THREE.Group,
        text: string,
        pos: THREE.Vector3,
        color: number
    ) {
        if (!this.font) return;
        const geo = new TextGeometry(text, {
            font: this.font,
            size: this.params.fontSize,
            height: 0.001,
        });
        const mat = new THREE.MeshBasicMaterial({ color });
        const mesh = new THREE.Mesh(geo, mat);
        geo.computeBoundingBox();
        if (geo.boundingBox) {
            const center = new THREE.Vector3();
            geo.boundingBox.getCenter(center);
            mesh.position.copy(pos.clone().sub(center.multiplyScalar(0.5)));
        } else {
            mesh.position.copy(pos);
        }
        group.add(mesh);
    }

    private niceStep(step: number): number {
        const pow10 = Math.pow(10, Math.floor(Math.log10(step)));
        const n = step / pow10;
        const nice = n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10;
        return nice * pow10;
    }

    private formatNumber(x: number): string {
        const abs = Math.abs(x);
        if (abs >= 1) return x.toFixed(1);
        if (abs >= 0.01) return x.toFixed(2);
        return x.toExponential(1);
    }

    public update() {
        if (!this.font) return;
        const box = new THREE.Box3().setFromObject(this.target);
        this.boxHelper.box.copy(box);
    }

    public dispose() {
        this.clear();
    }
}
