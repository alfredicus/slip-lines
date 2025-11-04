import * as THREE from 'three';

export interface ImageBackgroundConfig {
    imagePath?: string;
    gradientTop?: string;
    gradientBottom?: string;
    useGradient?: boolean;
}

/**
 * Creates a fixed background that doesn't rotate with the camera.
 * The background can be either an image or a generated gradient.
 */
export class FixedImageBackground {
    private scene: THREE.Scene;
    private texture: THREE.Texture | null = null;

    constructor(scene: THREE.Scene, config?: ImageBackgroundConfig) {
        this.scene = scene;

        if (config?.useGradient) {
            // Generate gradient as texture
            const gradientTexture = this.createGradientTexture(
                config.gradientTop || '#87CEEB',
                config.gradientBottom || '#E0F6FF'
            );
            this.applyBackground(gradientTexture);
        } else if (config?.imagePath) {
            // Load image from path
            this.loadImage(config.imagePath);
        }
    }

    /**
     * Creates a gradient texture using canvas
     */
    private createGradientTexture(
        topColor: string,
        bottomColor: string,
        width: number = 1024,
        height: number = 1024
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

    /**
     * Creates a radial gradient texture
     */
    createRadialGradient(
        centerColor: string,
        edgeColor: string,
        width: number = 1024,
        height: number = 1024
    ): THREE.CanvasTexture {
        const canvas = document.createElement('canvas');
        canvas.width = width;
        canvas.height = height;

        const context = canvas.getContext('2d')!;
        const centerX = width / 2;
        const centerY = height / 2;
        const radius = Math.min(width, height) / 2;

        const gradient = context.createRadialGradient(
            centerX, centerY, 0,
            centerX, centerY, radius
        );

        gradient.addColorStop(0, centerColor);
        gradient.addColorStop(1, edgeColor);

        context.fillStyle = gradient;
        context.fillRect(0, 0, width, height);

        const texture = new THREE.CanvasTexture(canvas);
        texture.needsUpdate = true;

        return texture;
    }

    /**
     * Creates a multi-stop gradient texture
     */
    createMultiStopGradient(
        colorStops: Array<{ offset: number; color: string }>,
        width: number = 1024,
        height: number = 1024,
        horizontal: boolean = false
    ): THREE.CanvasTexture {
        const canvas = document.createElement('canvas');
        canvas.width = width;
        canvas.height = height;

        const context = canvas.getContext('2d')!;
        const gradient = horizontal
            ? context.createLinearGradient(0, 0, width, 0)
            : context.createLinearGradient(0, 0, 0, height);

        colorStops.forEach(stop => {
            gradient.addColorStop(stop.offset, stop.color);
        });

        context.fillStyle = gradient;
        context.fillRect(0, 0, width, height);

        const texture = new THREE.CanvasTexture(canvas);
        texture.needsUpdate = true;

        return texture;
    }

    /**
     * Loads an image from a URL or path
     */
    loadImage(imagePath: string): Promise<void> {
        return new Promise((resolve, reject) => {
            const loader = new THREE.TextureLoader();
            loader.load(
                imagePath,
                (texture) => {
                    this.applyBackground(texture);
                    resolve();
                },
                undefined,
                (error) => {
                    console.error('Error loading background image:', error);
                    reject(error);
                }
            );
        });
    }

    /**
     * Applies a texture as the scene background
     * The texture will remain fixed and not rotate with the camera
     */
    private applyBackground(texture: THREE.Texture): void {
        this.texture = texture;
        this.scene.background = texture;
    }

    /**
     * Updates the background with a new gradient
     */
    updateGradient(topColor: string, bottomColor: string): void {
        const gradientTexture = this.createGradientTexture(topColor, bottomColor);
        this.applyBackground(gradientTexture);
    }

    /**
     * Updates the background with a new image
     */
    async updateImage(imagePath: string): Promise<void> {
        await this.loadImage(imagePath);
    }

    /**
     * Updates the background with a radial gradient
     */
    updateRadialGradient(centerColor: string, edgeColor: string): void {
        const texture = this.createRadialGradient(centerColor, edgeColor);
        this.applyBackground(texture);
    }

    /**
     * Applies a custom texture
     */
    applyCustomTexture(texture: THREE.Texture): void {
        this.applyBackground(texture);
    }

    /**
     * Removes the background
     */
    clear(): void {
        this.scene.background = null;
    }

    /**
     * Sets a solid color background
     */
    setSolidColor(color: string | number): void {
        this.scene.background = new THREE.Color(color);
    }

    /**
     * Cleanup
     */
    dispose(): void {
        if (this.texture) {
            this.texture.dispose();
            this.texture = null;
        }
        this.scene.background = null;
    }
}

/**
 * Preset gradient configurations
 */
export const BACKGROUND_PRESETS = {
    sky: { top: '#87CEEB', bottom: '#E0F6FF' },
    sunset: { top: '#FF6B35', bottom: '#F7931E' },
    ocean: { top: '#003973', bottom: '#0891B2' },
    forest: { top: '#1A4D2E', bottom: '#4E9F3D' },
    midnight: { top: '#0F0F1E', bottom: '#1E293B' },
    desert: { top: '#FFA500', bottom: '#FFD700' },
    arctic: { top: '#DBEAFE', bottom: '#F0F9FF' },
    volcanic: { top: '#DC2626', bottom: '#1C1917' },
    twilight: { top: '#4A148C', bottom: '#FF6F00' },
    dawn: { top: '#FFF9C4', bottom: '#FFAB91' }
};

/**
 * Helper function to create a complex gradient background
 */
export function createComplexGradient(
    scene: THREE.Scene,
    type: 'sunrise' | 'sunset' | 'northern-lights' | 'tropical' | 'grayscale'
): FixedImageBackground {
    const background = new FixedImageBackground(scene);

    switch (type) {
        case 'sunrise':
            background.applyCustomTexture(
                background.createMultiStopGradient([
                    { offset: 0.0, color: '#001F3F' },
                    { offset: 0.3, color: '#FF6B35' },
                    { offset: 0.6, color: '#FFA500' },
                    { offset: 1.0, color: '#FFE5B4' }
                ])
            );
            break;

        case 'sunset':
            background.applyCustomTexture(
                background.createMultiStopGradient([
                    { offset: 0.0, color: '#2C1B47' },
                    { offset: 0.4, color: '#D84315' },
                    { offset: 0.7, color: '#FF6F00' },
                    { offset: 1.0, color: '#FFB74D' }
                ])
            );
            break;

        case 'northern-lights':
            background.applyCustomTexture(
                background.createMultiStopGradient([
                    { offset: 0.0, color: '#0A1929' },
                    { offset: 0.3, color: '#1E3A5F' },
                    { offset: 0.5, color: '#00BCD4' },
                    { offset: 0.7, color: '#4CAF50' },
                    { offset: 1.0, color: '#0D1B2A' }
                ])
            );
            break;

        case 'tropical':
            background.applyCustomTexture(
                background.createMultiStopGradient([
                    { offset: 0.0, color: '#00ACC1' },
                    { offset: 0.5, color: '#26C6DA' },
                    { offset: 1.0, color: '#B2EBF2' }
                ])
            );
            break;

        case 'grayscale':
            background.applyCustomTexture(
                background.createMultiStopGradient([
                    { offset: 0.0, color: '#8d8c8cff' },
                    { offset: 1.0, color: '#000000' }
                ])
            );
            break;
    }

    return background;
}