import { ColorStop } from './ColorScale';
import { ColorMap, colorMapNames } from './colorMap';

/**
 * Converts a hex color number to a hex string
 */
function hexNumberToString(hexNumber: number): string {
    return '#' + hexNumber.toString(16).padStart(6, '0');
}

/**
 * Converts a ColorMap definition to ColorScale color stops
 * @param colorMapName - Name of the color map from colorMap.ts
 * @returns Array of ColorStop objects
 */
export function colorMapToColorStops(colorMapName: string): ColorStop[] {
    // Create a ColorMap instance with the specified name
    const colorMap = new ColorMap(colorMapName, 256);

    // Get the raw map colors
    const mapColors = colorMap.mapColors;

    // Convert to ColorStop format
    const colorStops: ColorStop[] = mapColors.map((entry: [number, number]) => ({
        position: entry[0],
        color: hexNumberToString(entry[1])
    }));

    return colorStops;
}

/**
 * Gets all available color map names
 */
export function getAvailableColorMaps(): string[] {
    return colorMapNames() as string[];
}

/**
 * Converts a ColorMap to a canvas gradient for preview
 */
export function createColorMapPreview(
    colorMapName: string,
    width: number = 256,
    height: number = 30
): HTMLCanvasElement {
    const colorMap = new ColorMap(colorMapName, 256);
    const canvas = colorMap.createCanvas(document, width);
    canvas.height = height;

    // Stretch the gradient to fill the height
    const ctx = canvas.getContext('2d')!;
    const imageData = ctx.getImageData(0, 0, width, 1);

    for (let y = 0; y < height; y++) {
        ctx.putImageData(imageData, 0, y);
    }

    return canvas;
}

/**
 * Gets a color from a ColorMap at a specific value
 */
export function getColorFromColorMap(
    colorMapName: string,
    value: number,
    min: number = 0,
    max: number = 1
): string {
    const colorMap = new ColorMap(colorMapName, 256);
    colorMap.setMin(min);
    colorMap.setMax(max);

    const color = colorMap.getColor(value);
    return `rgb(${Math.round(color.r * 255)}, ${Math.round(color.g * 255)}, ${Math.round(color.b * 255)})`;
}

/**
 * Interpolates between color stops at a given position
 */
export function interpolateColorStops(
    colorStops: ColorStop[],
    position: number
): string {
    const sortedStops = [...colorStops].sort((a, b) => a.position - b.position);

    if (position <= sortedStops[0].position) {
        return sortedStops[0].color;
    }

    if (position >= sortedStops[sortedStops.length - 1].position) {
        return sortedStops[sortedStops.length - 1].color;
    }

    for (let i = 0; i < sortedStops.length - 1; i++) {
        const stop1 = sortedStops[i];
        const stop2 = sortedStops[i + 1];

        if (position >= stop1.position && position <= stop2.position) {
            const t = (position - stop1.position) / (stop2.position - stop1.position);
            return interpolateHexColors(stop1.color, stop2.color, t);
        }
    }

    return sortedStops[0].color;
}

/**
 * Interpolates between two hex colors
 */
function interpolateHexColors(color1: string, color2: string, t: number): string {
    const c1 = hexToRgb(color1);
    const c2 = hexToRgb(color2);

    if (!c1 || !c2) return color1;

    const r = Math.round(c1.r + (c2.r - c1.r) * t);
    const g = Math.round(c1.g + (c2.g - c1.g) * t);
    const b = Math.round(c1.b + (c2.b - c1.b) * t);

    return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Converts hex color to RGB
 */
function hexToRgb(hex: string): { r: number, g: number, b: number } | null {
    // Remove # if present
    hex = hex.replace('#', '');

    const result = /^([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}

/**
 * Creates a custom color map with smooth interpolation
 */
export function createSmoothColorMap(
    colorMapName: string,
    numStops: number = 10
): ColorStop[] {
    const colorMap = new ColorMap(colorMapName, 256);
    const mapColors = colorMap.mapColors;

    const colorStops: ColorStop[] = [];

    for (let i = 0; i < numStops; i++) {
        const position = i / (numStops - 1);

        // Find which segment this position falls into
        let color = mapColors[0][1];

        for (let j = 0; j < mapColors.length - 1; j++) {
            if (position >= mapColors[j][0] && position <= mapColors[j + 1][0]) {
                const t = (position - mapColors[j][0]) / (mapColors[j + 1][0] - mapColors[j][0]);
                const c1 = hexNumberToString(mapColors[j][1]);
                const c2 = hexNumberToString(mapColors[j + 1][1]);
                color = interpolateHexColors(c1, c2, t);
                break;
            }
        }

        colorStops.push({
            position,
            color: typeof color === 'number' ? hexNumberToString(color) : color
        });
    }

    return colorStops;
}