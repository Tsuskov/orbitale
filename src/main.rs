use std::f64::consts::PI;

// Bohr radius in units of a₀ = 1
const BOHR_RADIUS: f64 = 1.0;

// Hydrogen wavefunction ψ_{n,l,m}(r, θ, φ)
fn hydrogen_wavefunction(r: f64, theta: f64, phi: f64, n: u32, l: u32, m: i32) -> f64 {
    if l >= n {
        return 0.0;
    }

    let rho = 2.0 * r / (n as f64 * BOHR_RADIUS);

    // Radial part: R_{n,l}(r)
    let radial = radial_part(rho, n, l);

    // Angular part: Y_{l,m}(θ, φ)
    let angular = angular_part(theta, phi, l, m);

    radial * angular
}

// Radial wavefunction
fn radial_part(rho: f64, n: u32, l: u32) -> f64 {
    if rho > 100.0 {
        return 0.0; // Avoid overflow
    }

    let exp_factor = (-rho / 2.0).exp();
    let norm = normalization_factor(n, l);

    match (n, l) {
        (1, 0) => norm * 2.0 * exp_factor,
        (2, 0) => norm * (2.0 - rho) * exp_factor,
        (2, 1) => norm * rho * exp_factor,
        (3, 0) => norm * (6.0 - 6.0 * rho + rho * rho) * exp_factor,
        (3, 1) => norm * (4.0 * rho - rho * rho) * exp_factor,
        (3, 2) => norm * rho * rho * exp_factor,
        (4, 0) => norm * (24.0 - 36.0 * rho + 12.0 * rho * rho - rho * rho * rho) * exp_factor,
        (4, 1) => norm * (20.0 * rho - 10.0 * rho * rho + rho * rho * rho) * exp_factor,
        (4, 2) => norm * (6.0 * rho * rho - rho * rho * rho) * exp_factor,
        (4, 3) => norm * rho * rho * rho * exp_factor,
        _ => 0.0,
    }
}

// Normalization constants (simplified)
fn normalization_factor(n: u32, l: u32) -> f64 {
    match (n, l) {
        (1, 0) => 2.0,
        (2, 0) => 1.0 / (2.0 * 2.0_f64.sqrt()),
        (2, 1) => 1.0 / (2.0 * 6.0_f64.sqrt()),
        (3, 0) => 2.0 / (3.0 * 3.0_f64.sqrt()),
        (3, 1) => 4.0 / (27.0_f64.sqrt()),
        (3, 2) => 4.0 / (81.0 * 30.0_f64.sqrt()),
        (4, 0) => 1.0 / (4.0 * 6.0_f64.sqrt()),
        (4, 1) => 1.0 / (16.0 * 30.0_f64.sqrt()),
        (4, 2) => 1.0 / (64.0 * 35.0_f64.sqrt()),
        (4, 3) => 1.0 / (128.0 * 140.0_f64.sqrt()),
        _ => 1.0,
    }
}

// Simplified angular part (spherical harmonics)
fn angular_part(theta: f64, _phi: f64, l: u32, m: i32) -> f64 {
    match (l, m) {
        (0, 0) => 1.0 / (2.0 * PI.sqrt()),                 // s orbital
        (1, 0) => (3.0 / (4.0 * PI)).sqrt() * theta.cos(), // p_z
        (1, 1) | (1, -1) => (3.0 / (8.0 * PI)).sqrt() * theta.sin(), // p_x, p_y
        (2, 0) => (5.0 / (16.0 * PI)).sqrt() * (3.0 * theta.cos().powi(2) - 1.0), // d_z²
        (2, 1) | (2, -1) => (15.0 / (8.0 * PI)).sqrt() * theta.sin() * theta.cos(), // d_xz, d_yz
        (2, 2) | (2, -2) => (15.0 / (32.0 * PI)).sqrt() * theta.sin().powi(2), // d_x²-y², d_xy
        (3, 0) => (7.0 / (16.0 * PI)).sqrt() * (5.0 * theta.cos().powi(3) - 3.0 * theta.cos()), // f_z³
        (3, 1) | (3, -1) => {
            (21.0 / (64.0 * PI)).sqrt() * theta.sin() * (5.0 * theta.cos().powi(2) - 1.0)
        } // f_xz², f_yz²
        (3, 2) | (3, -2) => (105.0 / (32.0 * PI)).sqrt() * theta.sin().powi(2) * theta.cos(), // f_xyz, f_...
        (3, 3) | (3, -3) => (35.0 / (64.0 * PI)).sqrt() * theta.sin().powi(3), // f_x³, f_y³
        _ => 0.0,
    }
}

// Electron density
fn electron_density(r: f64, theta: f64, phi: f64, n: u32, l: u32, m: i32) -> f64 {
    let psi = hydrogen_wavefunction(r, theta, phi, n, l, m);
    (psi * psi).abs()
}

// Render 2D slice in xy-plane (z=0, which is θ=π/2)
fn render_orbital(n: u32, l: u32, m: i32, label: &str) {
    println!("\n=== {} orbital (n={}, l={}, m={}) ===", label, n, l, m);

    let width = 60;
    let height = 30;
    let max_r = (n as f64) * (n as f64) * 3.0; // Scale with n²

    let mut density_map: Vec<Vec<f64>> = vec![vec![0.0; width]; height]; // density_map[row][col]

    // Calculate unnormalized density values throuout map
    let mut max_density = 0.0;
    for row in 0..height {
        for col in 0..width {
            let x = (col as f64 - width as f64 / 2.0) * max_r / (width as f64 / 2.0);
            let y = (height as f64 / 2.0 - row as f64) * max_r / (height as f64 / 2.0);
            let r = (x * x + y * y).sqrt();
            let theta = if r > 1e-6 {
                y.atan2(x) + PI / 2.0
            } else {
                PI / 2.0
            };

            let density = electron_density(r, theta, 0.0, n, l, m);

            if density > max_density {
                max_density = density;
            } // determine max density for later normalization

            density_map[row][col] = density;
        }
    }

    let scale = if max_density > 0.0 {
        30.0 / max_density
    } else {
        1.0
    };

    for row in 0..height {
        for col in 0..width {
            let scaled = (density_map[row][col] * scale) as u32; // normalized density value

            let char = match scaled {
                0..=1 => ' ',
                2..=5 => '·',
                6..=10 => '∘',
                11..=15 => '◯',
                16..=22 => '●',
                _ => '█',
            };

            print!("{}", char);
        }
        println!();
    }
}

fn main() {
    println!("Hydrogen Orbital Visualizer (2D xy-plane cross-section)");
    println!("========================================================\n");

    // s orbitals
    println!("\n--- S ORBITALS ---");
    render_orbital(1, 0, 0, "1s");
    render_orbital(2, 0, 0, "2s");
    render_orbital(3, 0, 0, "3s");

    // p orbitals
    println!("\n--- P ORBITALS ---");
    render_orbital(2, 1, 0, "2p");
    render_orbital(3, 1, 0, "3p");
    render_orbital(4, 1, 0, "4p");

    // d orbitals
    println!("\n--- D ORBITALS ---");
    render_orbital(3, 2, 0, "3d");
    render_orbital(4, 2, 0, "4d");

    // f orbitals
    println!("\n--- F ORBITALS ---");
    render_orbital(4, 3, 0, "4f");

    println!("\nDensity scale: · ∘ ◯ ● █ (increasing)");
}
