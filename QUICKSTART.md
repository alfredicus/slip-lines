# Quick Start Guide - 3D Slip Lines Visualization

## 🚀 Get Started in 3 Minutes

### Step 1: Install Dependencies
```bash
npm install
```

### Step 2: Run Development Server
```bash
npm run dev
```

Your browser will automatically open to `http://localhost:5173`

### Step 3: Start Exploring!

Use the GUI on the right side to:
- Adjust stress tensor parameters (σ_x, σ_y, σ_z)
- Toggle equipotentials, slip lines, and stress axes
- Change colors for each visualization element

## 📁 File Organization

```
project/
├── SlipLinesVisualization.tsx     ← Main component (TypeScript)
├── SlipLinesVisualization.jsx     ← Alternative (JavaScript)
├── main.tsx                        ← Entry point
├── index.html                      ← HTML template
├── index.css                       ← Global styles
├── vite.config.ts                 ← Build configuration
├── tsconfig.json                  ← TypeScript config
├── package.json                   ← Dependencies
└── README.md                       ← Full documentation
```

## 🎮 Mouse Controls

| Action | Effect |
|--------|--------|
| **Drag** | Rotate view |
| **Scroll** | Zoom in/out |

## 🎨 GUI Controls

### Stress Tensor Parameters
Three sliders to adjust the stress components:
- **σ_x**: Range from -5 to +5 (default: 1.0)
- **σ_y**: Range from -5 to +5 (default: 3.0)
- **σ_z**: Range from -5 to +5 (default: 2.0)

### Visualization Toggles
- **☐ Show Equipotentials** - Display constant stress curves
- **☐ Show Slip Lines** - Display theoretical slip lines
- **☐ Show Stress Axes** - Display principal stress directions (RGB = σ₁σ₂σ₃)

### Color Controls
- Click the color box next to each option to customize colors

## 💡 Example Configurations

### Pure Compression (Vertical)
```
σ_x = 3.0
σ_y = 1.0
σ_z = 1.0
```
Results in vertical slip lines and symmetric equipotentials

### Triaxial Extension
```
σ_x = -3.0
σ_y = 1.0
σ_z = 0.5
```
Results in extensional slip patterns

### Strike-Slip Regime
```
σ_x = 2.0
σ_y = 3.0
σ_z = 1.5
```
Results in lateral slip lines

## 🔧 Customization

### Change Default Values
Edit the initial state in `SlipLinesVisualization.tsx`:

```tsx
const [controls, setControls] = useState<ControlState>({
  sigmaX: 1.0,      // ← Change here
  sigmaY: 3.0,      // ← Change here
  sigmaZ: 2.0,      // ← Change here
  // ...
});
```

### Adjust Curve Density
Look for these variables and modify:

```tsx
const nCurves = 12;           // Slip lines per quadrant (more = denser)
const nIsoCurves = 8;         // Equipotential levels (more = more detail)
const curveResolution = 200;  // Points per curve (more = smoother)
```

### Change Sphere Appearance
Find the sphere material creation and modify:

```tsx
const sphereMaterial = new THREE.MeshPhongMaterial({
  color: 0xe8e8e8,        // Change hex color
  opacity: 0.3,           // 0 = transparent, 1 = opaque
  wireframe: false        // Set true to see mesh
});
```

## 🚢 Build for Production

```bash
npm run build
```

Output files will be in the `dist/` folder. Deploy this folder to any static hosting service:
- GitHub Pages
- Netlify
- Vercel
- AWS S3
- etc.

## 🐛 Troubleshooting

### Port 5173 Already in Use?
```bash
npm run dev -- --port 3000
```

### Component Not Rendering?
1. Check browser console for errors (F12)
2. Ensure WebGL is enabled in your browser
3. Try a different browser

### Slow Performance?
1. Reduce `curveResolution` from 200 to 100
2. Reduce `nCurves` from 12 to 6
3. Disable stress axes visualization

### Colors Not Changing?
1. Make sure to use valid hex color codes (#RRGGBB)
2. Try closing and reopening the browser

## 📚 Learn More

### Three.js Documentation
- Official site: https://threejs.org
- Examples: https://threejs.org/examples
- Documentation: https://threejs.org/docs

### React Hooks
- https://react.dev/reference/react

### TypeScript
- https://www.typescriptlang.org

## 🤝 Extending the Project

### Add DEM Data Visualization
1. Implement `loadDEMData()` function
2. Generate trajectory curves
3. Add visualization toggle in GUI

### Add Export Functionality
```tsx
const exportImage = () => {
  const link = document.createElement('a');
  link.href = renderer.domElement.toDataURL();
  link.download = 'slip-lines.png';
  link.click();
};
```

### Add Animation
```tsx
useEffect(() => {
  let angle = 0;
  const interval = setInterval(() => {
    angle += 0.01;
    camera.position.applyAxisAngle(new THREE.Vector3(0, 1, 0), 0.01);
    camera.lookAt(0, 0, 0);
  }, 50);
  return () => clearInterval(interval);
}, []);
```

## 📦 Dependencies

| Package | Purpose |
|---------|---------|
| `react` | UI framework |
| `react-dom` | DOM rendering |
| `three` | 3D graphics |
| `typescript` | Type safety |
| `vite` | Build tool |

## ✨ Features Implemented

✅ 3D sphere with semi-transparent material
✅ Equipotential curves generation
✅ Slip lines generation
✅ Principal stress axes visualization
✅ Interactive mouse controls (rotate, zoom)
✅ Real-time parameter adjustment
✅ Color customization
✅ TypeScript support
✅ Responsive layout
✅ Smooth rendering at 60fps

## 📝 License

This project is a TypeScript/React conversion of the original Python implementation.

---

**Need help?** Check the full README.md for more detailed information!
