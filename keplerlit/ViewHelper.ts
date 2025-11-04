import {
    CylinderGeometry,
    CanvasTexture,
    Color,
    Euler,
    Mesh,
    MeshBasicMaterial,
    Object3D,
    OrthographicCamera,
    Quaternion,
    Raycaster,
    Sprite,
    SpriteMaterial,
    SRGBColorSpace,
    Vector2,
    Vector3,
    Vector4,
    Camera,
    WebGLRenderer
} from 'three';

interface LabelOptions {
    labelX?: string;
    labelY?: string;
    labelZ?: string;
    font?: string;
    color?: string;
    radius?: number;
}

type AxisType = 'posX' | 'posY' | 'posZ' | 'negX' | 'negY' | 'negZ';

interface AxisUserData {
    type: AxisType;
}

/**
 * A special type of helper that visualizes the camera's transformation
 * in a small viewport area as an axes helper. Such a helper is often wanted
 * in 3D modeling tools or scene editors like the [three.js editor]{@link https://threejs.org/editor}.
 *
 * The helper allows to click on the X, Y and Z axes which animates the camera
 * so it looks along the selected axis.
 *
 * @augments Object3D
 * @three_import import { ViewHelper } from 'three/addons/helpers/ViewHelper.js';
 */
class ViewHelper extends Object3D {
    /**
     * This flag can be used for type testing.
     *
     * @type {boolean}
     * @readonly
     * @default true
     */
    readonly isViewHelper: boolean = true;

    /**
     * Whether the helper is currently animating or not.
     *
     * @type {boolean}
     * @readonly
     * @default false
     */
    animating: boolean = false;

    /**
     * The helper's center point.
     *
     * @type {Vector3}
     */
    center: Vector3;

    private camera: Camera;
    private domElement: HTMLElement;
    private options: LabelOptions;
    private interactiveObjects: Sprite[];
    private raycaster: Raycaster;
    private mouse: Vector2;
    private dummy: Object3D;
    private orthoCamera: OrthographicCamera;
    private geometry: CylinderGeometry;
    private xAxis: Mesh;
    private yAxis: Mesh;
    private zAxis: Mesh;
    private posXAxisHelper: Sprite;
    private posYAxisHelper: Sprite;
    private posZAxisHelper: Sprite;
    private negXAxisHelper: Sprite;
    private negYAxisHelper: Sprite;
    private negZAxisHelper: Sprite;
    private targetPosition: Vector3;
    private targetQuaternion: Quaternion;
    private q1: Quaternion;
    private q2: Quaternion;
    private viewport: Vector4;
    private radius: number;
    private point: Vector3;
    private dim: number;
    private turnRate: number;
    private color1: Color;
    private color2: Color;
    private color3: Color;
    private color4: Color;

    /**
     * Constructs a new view helper.
     *
     * @param {Camera} camera - The camera whose transformation should be visualized.
     * @param {HTMLElement} domElement - The DOM element that is used to render the view.
     */
    constructor(camera: Camera, domElement: HTMLElement) {
        super();

        this.camera = camera;
        this.domElement = domElement;
        this.center = new Vector3();

        this.color1 = new Color('#ff4466');
        this.color2 = new Color('#88ff44');
        this.color3 = new Color('#4488ff');
        this.color4 = new Color('#000000');

        this.options = {};
        this.interactiveObjects = [];
        this.raycaster = new Raycaster();
        this.mouse = new Vector2();
        this.dummy = new Object3D();

        this.orthoCamera = new OrthographicCamera(-2, 2, 2, -2, 0, 4);
        this.orthoCamera.position.set(0, 0, 2);

        this.geometry = new CylinderGeometry(0.04, 0.04, 0.8, 5)
            .rotateZ(-Math.PI / 2)
            .translate(0.4, 0, 0);

        this.xAxis = new Mesh(this.geometry, this.getAxisMaterial(this.color1));
        this.yAxis = new Mesh(this.geometry, this.getAxisMaterial(this.color2));
        this.zAxis = new Mesh(this.geometry, this.getAxisMaterial(this.color3));

        this.yAxis.rotation.z = Math.PI / 2;
        this.zAxis.rotation.y = -Math.PI / 2;

        this.add(this.xAxis);
        this.add(this.zAxis);
        this.add(this.yAxis);

        const spriteMaterial1 = this.getSpriteMaterial(this.color1);
        const spriteMaterial2 = this.getSpriteMaterial(this.color2);
        const spriteMaterial3 = this.getSpriteMaterial(this.color3);
        const spriteMaterial4 = this.getSpriteMaterial(this.color4);

        this.posXAxisHelper = new Sprite(spriteMaterial1);
        this.posYAxisHelper = new Sprite(spriteMaterial2);
        this.posZAxisHelper = new Sprite(spriteMaterial3);
        this.negXAxisHelper = new Sprite(spriteMaterial4);
        this.negYAxisHelper = new Sprite(spriteMaterial4);
        this.negZAxisHelper = new Sprite(spriteMaterial4);

        this.posXAxisHelper.position.x = 1;
        this.posYAxisHelper.position.y = 1;
        this.posZAxisHelper.position.z = 1;
        this.negXAxisHelper.position.x = -1;
        this.negYAxisHelper.position.y = -1;
        this.negZAxisHelper.position.z = -1;

        this.negXAxisHelper.material.opacity = 0.2;
        this.negYAxisHelper.material.opacity = 0.2;
        this.negZAxisHelper.material.opacity = 0.2;

        (this.posXAxisHelper.userData as AxisUserData).type = 'posX';
        (this.posYAxisHelper.userData as AxisUserData).type = 'posY';
        (this.posZAxisHelper.userData as AxisUserData).type = 'posZ';
        (this.negXAxisHelper.userData as AxisUserData).type = 'negX';
        (this.negYAxisHelper.userData as AxisUserData).type = 'negY';
        (this.negZAxisHelper.userData as AxisUserData).type = 'negZ';

        this.add(this.posXAxisHelper);
        this.add(this.posYAxisHelper);
        this.add(this.posZAxisHelper);
        this.add(this.negXAxisHelper);
        this.add(this.negYAxisHelper);
        this.add(this.negZAxisHelper);

        this.interactiveObjects.push(this.posXAxisHelper);
        this.interactiveObjects.push(this.posYAxisHelper);
        this.interactiveObjects.push(this.posZAxisHelper);
        this.interactiveObjects.push(this.negXAxisHelper);
        this.interactiveObjects.push(this.negYAxisHelper);
        this.interactiveObjects.push(this.negZAxisHelper);

        this.point = new Vector3();
        this.dim = 128;
        this.turnRate = 2 * Math.PI; // turn rate in angles per second

        this.targetPosition = new Vector3();
        this.targetQuaternion = new Quaternion();
        this.q1 = new Quaternion();
        this.q2 = new Quaternion();
        this.viewport = new Vector4();
        this.radius = 0;
    }

    /**
     * Renders the helper in a separate view in the bottom-right corner
     * of the viewport.
     *
     * @param {WebGLRenderer} renderer - The renderer.
     */
    render(renderer: WebGLRenderer): void {
        this.quaternion.copy(this.camera.quaternion).invert();
        this.updateMatrixWorld();

        this.point.set(0, 0, 1);
        this.point.applyQuaternion(this.camera.quaternion);

        const x = this.domElement.offsetWidth - this.dim;
        const y = (renderer as any).isWebGPURenderer
            ? this.domElement.offsetHeight - this.dim
            : 0;

        renderer.clearDepth();

        renderer.getViewport(this.viewport);
        renderer.setViewport(x, y, this.dim, this.dim);

        renderer.render(this, this.orthoCamera);

        renderer.setViewport(
            this.viewport.x,
            this.viewport.y,
            this.viewport.z,
            this.viewport.w
        );
    }

    /**
     * This method should be called when a click or pointer event
     * has happened in the app.
     *
     * @param {PointerEvent} event - The event to process.
     * @return {boolean} Whether an intersection with the helper has been detected or not.
     */
    handleClick(event: PointerEvent): boolean {
        if (this.animating === true) return false;

        const rect = this.domElement.getBoundingClientRect();
        const offsetX = rect.left + (this.domElement.offsetWidth - this.dim);
        const offsetY = rect.top + (this.domElement.offsetHeight - this.dim);
        this.mouse.x = ((event.clientX - offsetX) / (rect.right - offsetX)) * 2 - 1;
        this.mouse.y = -((event.clientY - offsetY) / (rect.bottom - offsetY)) * 2 + 1;

        this.raycaster.setFromCamera(this.mouse, this.orthoCamera);

        const intersects = this.raycaster.intersectObjects(this.interactiveObjects);

        if (intersects.length > 0) {
            const intersection = intersects[0];
            const object = intersection.object as Sprite;

            this.prepareAnimationData(object, this.center);

            this.animating = true;

            return true;
        } else {
            return false;
        }
    }

    /**
     * Sets labels for each axis. By default, they are unlabeled.
     *
     * @param {string|undefined} labelX - The label for the x-axis.
     * @param {string|undefined} labelY - The label for the y-axis.
     * @param {string|undefined} labelZ - The label for the z-axis.
     */
    setLabels(labelX?: string, labelY?: string, labelZ?: string): void {
        this.options.labelX = labelX;
        this.options.labelY = labelY;
        this.options.labelZ = labelZ;

        this.updateLabels();
    }

    /**
     * Sets the label style. Has no effect when the axes are unlabeled.
     *
     * @param {string} [font='24px Arial'] - The label font.
     * @param {string} [color='#000000'] - The label color.
     * @param {number} [radius=14] - The label radius.
     */
    setLabelStyle(font?: string, color?: string, radius?: number): void {
        this.options.font = font;
        this.options.color = color;
        this.options.radius = radius;

        this.updateLabels();
    }

    /**
     * Updates the helper. This method should be called in the app's animation
     * loop.
     *
     * @param {number} delta - The delta time in seconds.
     */
    update(delta: number): void {
        const step = delta * this.turnRate;

        // animate position by doing a slerp and then scaling the position on the unit sphere
        this.q1.rotateTowards(this.q2, step);
        this.camera.position
            .set(0, 0, 1)
            .applyQuaternion(this.q1)
            .multiplyScalar(this.radius)
            .add(this.center);

        // animate orientation
        this.camera.quaternion.rotateTowards(this.targetQuaternion, step);

        if (this.q1.angleTo(this.q2) === 0) {
            this.animating = false;
        }
    }

    /**
     * Frees the GPU-related resources allocated by this instance. Call this
     * method whenever this instance is no longer used in your app.
     */
    dispose(): void {
        this.geometry.dispose();

        // this.xAxis.material.dispose();
        // this.yAxis.material.dispose();
        // this.zAxis.material.dispose();

        this.posXAxisHelper.material.map?.dispose();
        this.posYAxisHelper.material.map?.dispose();
        this.posZAxisHelper.material.map?.dispose();
        this.negXAxisHelper.material.map?.dispose();
        this.negYAxisHelper.material.map?.dispose();
        this.negZAxisHelper.material.map?.dispose();

        this.posXAxisHelper.material.dispose();
        this.posYAxisHelper.material.dispose();
        this.posZAxisHelper.material.dispose();
        this.negXAxisHelper.material.dispose();
        this.negYAxisHelper.material.dispose();
        this.negZAxisHelper.material.dispose();
    }

    private prepareAnimationData(object: Sprite, focusPoint: Vector3): void {
        const userData = object.userData as AxisUserData;

        switch (userData.type) {
            case 'posX':
                this.targetPosition.set(1, 0, 0);
                this.targetQuaternion.setFromEuler(new Euler(0, Math.PI * 0.5, 0));
                break;

            case 'posY':
                this.targetPosition.set(0, 1, 0);
                this.targetQuaternion.setFromEuler(new Euler(-Math.PI * 0.5, 0, 0));
                break;

            case 'posZ':
                this.targetPosition.set(0, 0, 1);
                this.targetQuaternion.setFromEuler(new Euler());
                break;

            case 'negX':
                this.targetPosition.set(-1, 0, 0);
                this.targetQuaternion.setFromEuler(new Euler(0, -Math.PI * 0.5, 0));
                break;

            case 'negY':
                this.targetPosition.set(0, -1, 0);
                this.targetQuaternion.setFromEuler(new Euler(Math.PI * 0.5, 0, 0));
                break;

            case 'negZ':
                this.targetPosition.set(0, 0, -1);
                this.targetQuaternion.setFromEuler(new Euler(0, Math.PI, 0));
                break;

            default:
                console.error('ViewHelper: Invalid axis.');
        }

        this.radius = this.camera.position.distanceTo(focusPoint);
        this.targetPosition.multiplyScalar(this.radius).add(focusPoint);

        this.dummy.position.copy(focusPoint);

        this.dummy.lookAt(this.camera.position);
        this.q1.copy(this.dummy.quaternion);

        this.dummy.lookAt(this.targetPosition);
        this.q2.copy(this.dummy.quaternion);
    }

    private getAxisMaterial(color: Color): MeshBasicMaterial {
        return new MeshBasicMaterial({ color: color, toneMapped: false });
    }

    private getSpriteMaterial(color: Color, text?: string): SpriteMaterial {
        const { font = '24px Arial', color: labelColor = '#000000', radius = 14 } = this.options;

        const canvas = document.createElement('canvas');
        canvas.width = 64;
        canvas.height = 64;

        const context = canvas.getContext('2d')!;
        context.beginPath();
        context.arc(32, 32, radius, 0, 2 * Math.PI);
        context.closePath();
        context.fillStyle = color.getStyle();
        context.fill();

        if (text) {
            context.font = font;
            context.textAlign = 'center';
            context.fillStyle = labelColor;
            context.fillText(text, 32, 41);
        }

        const texture = new CanvasTexture(canvas);
        texture.colorSpace = SRGBColorSpace;

        return new SpriteMaterial({ map: texture, toneMapped: false });
    }

    private updateLabels(): void {
        this.posXAxisHelper.material.map?.dispose();
        this.posYAxisHelper.material.map?.dispose();
        this.posZAxisHelper.material.map?.dispose();

        this.posXAxisHelper.material.dispose();
        this.posYAxisHelper.material.dispose();
        this.posZAxisHelper.material.dispose();

        this.posXAxisHelper.material = this.getSpriteMaterial(this.color1, this.options.labelX);
        this.posYAxisHelper.material = this.getSpriteMaterial(this.color2, this.options.labelY);
        this.posZAxisHelper.material = this.getSpriteMaterial(this.color3, this.options.labelZ);
    }
}

export { ViewHelper };