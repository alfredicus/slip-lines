import * as THREE from 'three';

export interface GradientConfig {
    topColor: string;
    bottomColor: string;
    offset?: number;
    exponent?: number;
}

export type GradientPreset =
    | 'sky'
    | 'sunset'
    | 'ocean'
    | 'forest'
    | 'midnight'
    | 'desert'
    | 'arctic'
    | 'volcanic';

export const GRADIENT_PRESETS: Record<GradientPreset, GradientConfig> = {
    sky: {
        topColor: '#87CEEB',
        bottomColor: '#E0F6FF',
        offset: 0.5,
        exponent: 0.6
    },
    sunset: {
        topColor: '#FF6B35',
        bottomColor: '#F7931E',
        offset: 0.4,
        exponent: 1.2
    },
    ocean: {
        topColor: '#003973',
        bottomColor: '#0891B2',
        offset: 0.5,
        exponent: 0.8
    },
    forest: {
        topColor: '#1A4D2E',
        bottomColor: '#4E9F3D',
        offset: 0.5,
        exponent: 1.0
    },
    midnight: {
        topColor: '#0F0F1E',
        bottomColor: '#1E293B',
        offset: 0.5,
        exponent: 0.7
    },
    desert: {
        topColor: '#FFA500',
        bottomColor: '#FFD700',
        offset: 0.45,
        exponent: 1.5
    },
    arctic: {
        topColor: '#DBEAFE',
        bottomColor: '#F0F9FF',
        offset: 0.5,
        exponent: 0.5
    },
    volcanic: {
        topColor: '#DC2626',
        bottomColor: '#1C1917',
        offset: 0.3,
        exponent: 2.0
    }
};

export class GradientBackground {
    private scene: THREE.Scene;
    private mesh: THREE.Mesh | null = null;
    private config: Required<GradientConfig>;

    constructor(scene: THREE.Scene, config: GradientConfig | GradientPreset = 'sky') {
        this.scene = scene;

        // Handle preset or custom config
        const gradientConfig = typeof config === 'string'
            ? GRADIENT_PRESETS[config]
            : config;

        this.config = {
            topColor: gradientConfig.topColor,
            bottomColor: gradientConfig.bottomColor,
            offset: gradientConfig.offset ?? 0.5,
            exponent: gradientConfig.exponent ?? 1.0
        };

        this.create();
    }

    /**
     * Creates the gradient background using a shader material
     */
    private create(): void {
        const vertexShader = `
			varying vec3 vWorldPosition;

			void main() {
				vec4 worldPosition = modelMatrix * vec4(position, 1.0);
				vWorldPosition = worldPosition.xyz;
				gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
			}
		`;

        const fragmentShader = `
			uniform vec3 topColor;
			uniform vec3 bottomColor;
			uniform float offset;
			uniform float exponent;

			varying vec3 vWorldPosition;

			void main() {
				float h = normalize(vWorldPosition + offset).y;
				float gradient = clamp(pow(max(h, 0.0), exponent), 0.0, 1.0);
				gl_FragColor = vec4(mix(bottomColor, topColor, gradient), 1.0);
			}
		`;

        const uniforms = {
            topColor: { value: new THREE.Color(this.config.topColor) },
            bottomColor: { value: new THREE.Color(this.config.bottomColor) },
            offset: { value: this.config.offset },
            exponent: { value: this.config.exponent }
        };

        const skyGeo = new THREE.SphereGeometry(500, 32, 15);
        const skyMat = new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: vertexShader,
            fragmentShader: fragmentShader,
            side: THREE.BackSide
        });

        this.mesh = new THREE.Mesh(skyGeo, skyMat);
        this.scene.add(this.mesh);
    }

    /**
     * Updates the gradient colors
     */
    updateColors(topColor: string, bottomColor: string): void {
        if (!this.mesh) return;

        this.config.topColor = topColor;
        this.config.bottomColor = bottomColor;

        const material = this.mesh.material as THREE.ShaderMaterial;
        material.uniforms.topColor.value.set(topColor);
        material.uniforms.bottomColor.value.set(bottomColor);
    }

    /**
     * Updates the gradient offset (vertical position)
     */
    updateOffset(offset: number): void {
        if (!this.mesh) return;

        this.config.offset = offset;
        const material = this.mesh.material as THREE.ShaderMaterial;
        material.uniforms.offset.value = offset;
    }

    /**
     * Updates the gradient exponent (smoothness/sharpness)
     */
    updateExponent(exponent: number): void {
        if (!this.mesh) return;

        this.config.exponent = exponent;
        const material = this.mesh.material as THREE.ShaderMaterial;
        material.uniforms.exponent.value = exponent;
    }

    /**
     * Applies a preset gradient
     */
    applyPreset(preset: GradientPreset): void {
        const config = GRADIENT_PRESETS[preset];
        this.updateColors(config.topColor, config.bottomColor);
        this.updateOffset(config.offset ?? 0.5);
        this.updateExponent(config.exponent ?? 1.0);
    }

    /**
     * Gets the current configuration
     */
    getConfig(): Required<GradientConfig> {
        return { ...this.config };
    }

    /**
     * Removes the gradient background from the scene
     */
    dispose(): void {
        if (!this.mesh) return;

        const material = this.mesh.material as THREE.ShaderMaterial;
        material.dispose();
        this.mesh.geometry.dispose();
        this.scene.remove(this.mesh);
        this.mesh = null;
    }
}

/**
 * Alternative: Simple CSS-based gradient background
 * This is the simplest approach and doesn't use Three.js rendering
 */
export function applyCSSGradient(
    container: HTMLElement,
    topColor: string,
    bottomColor: string
): void {
    container.style.background = `linear-gradient(to bottom, ${topColor}, ${bottomColor})`;
}

/**
 * Creates a canvas texture gradient for use as scene.background
 */
export function createCanvasGradient(
    topColor: string,
    bottomColor: string,
    width: number = 512,
    height: number = 512
): THREE.CanvasTexture {
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;

    const context = canvas.getContext('2d')!;
    const gradient = context.createLinearGradient(0, 0, 0, height);
    gradient.addColorStop(0, topColor);
    gradient.addColorStop(1, bottomColor);

    context.fillStyle = gradient;
    context.fillRect(0, 0, width, height);

    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;

    return texture;
}