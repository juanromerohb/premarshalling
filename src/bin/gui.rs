#![cfg_attr(target_os = "windows", windows_subsystem = "windows")]

//! Desktop GUI entrypoint.

use lct_rust::gui::app::CpmpGuiApp;

fn main() -> eframe::Result<()> {
    configure_linux_backend();

    let native_options = eframe::NativeOptions {
        renderer: eframe::Renderer::Glow,
        viewport: eframe::egui::ViewportBuilder::default()
            .with_title("CPMP-LCT Desktop GUI")
            .with_inner_size([1280.0, 900.0])
            .with_min_inner_size([980.0, 700.0]),
        ..Default::default()
    };

    if let Err(err) = eframe::run_native(
        "CPMP-LCT Desktop GUI",
        native_options,
        Box::new(|cc| Box::new(CpmpGuiApp::new(cc))),
    ) {
        eprintln!("GUI startup failed: {err}");
        eprintln!("Try running with:");
        eprintln!("  WINIT_UNIX_BACKEND=x11 LIBGL_ALWAYS_SOFTWARE=1 cargo run --bin gui");
        return Err(err);
    }

    Ok(())
}

fn configure_linux_backend() {
    #[cfg(target_os = "linux")]
    {
        if std::env::var_os("WINIT_UNIX_BACKEND").is_none() {
            std::env::set_var("WINIT_UNIX_BACKEND", "x11");
        }
        if std::env::var_os("LIBGL_ALWAYS_SOFTWARE").is_none() {
            std::env::set_var("LIBGL_ALWAYS_SOFTWARE", "1");
        }
    }
}
