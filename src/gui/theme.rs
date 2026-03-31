//! Shared visual theme helpers for the desktop GUI.

use eframe::egui;
use egui::{
    Color32, Context, FontFamily, FontId, Frame, Margin, Rounding, Stroke, TextStyle, Visuals,
};

/// Main window background color.
pub const BG: Color32 = Color32::from_rgb(16, 20, 28);
/// Primary panel background color.
pub const PANEL: Color32 = Color32::from_rgb(23, 30, 43);
/// Alternate panel/background accent color.
pub const PANEL_ALT: Color32 = Color32::from_rgb(30, 39, 55);
/// Primary accent color.
pub const ACCENT: Color32 = Color32::from_rgb(229, 138, 34);
/// Success/status-positive color.
pub const SUCCESS: Color32 = Color32::from_rgb(58, 179, 120);
/// Warning/status-pending color.
pub const WARNING: Color32 = Color32::from_rgb(224, 173, 63);
/// Error/status-negative color.
pub const DANGER: Color32 = Color32::from_rgb(211, 96, 96);
/// Border and separator color.
pub const BORDER: Color32 = Color32::from_rgb(66, 79, 100);

/// Applies the crate's custom theme to an egui context.
pub fn apply(ctx: &Context) {
    let mut visuals = Visuals::dark();
    visuals.window_fill = BG;
    visuals.panel_fill = BG;
    visuals.extreme_bg_color = PANEL;
    visuals.widgets.noninteractive.bg_fill = PANEL;
    visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, BORDER);
    visuals.widgets.inactive.bg_fill = PANEL_ALT;
    visuals.widgets.inactive.fg_stroke.color = Color32::from_rgb(230, 236, 246);
    visuals.widgets.active.bg_fill = Color32::from_rgb(60, 83, 114);
    visuals.widgets.hovered.bg_fill = Color32::from_rgb(74, 102, 140);
    visuals.selection.bg_fill = ACCENT;
    visuals.selection.stroke = Stroke::new(1.0, Color32::BLACK);

    let mut style = (*ctx.style()).clone();
    style.visuals = visuals;
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 6.0);
    style.text_styles = [
        (
            TextStyle::Heading,
            FontId::new(26.0, FontFamily::Proportional),
        ),
        (
            TextStyle::Name("Heading2".into()),
            FontId::new(22.0, FontFamily::Proportional),
        ),
        (TextStyle::Body, FontId::new(18.0, FontFamily::Proportional)),
        (
            TextStyle::Button,
            FontId::new(17.0, FontFamily::Proportional),
        ),
        (
            TextStyle::Monospace,
            FontId::new(16.0, FontFamily::Monospace),
        ),
        (
            TextStyle::Small,
            FontId::new(15.0, FontFamily::Proportional),
        ),
    ]
    .into();
    ctx.set_style(style);
    ctx.set_pixels_per_point(1.15);
}

/// Returns the standard card frame used across the GUI.
pub fn card_frame() -> Frame {
    Frame::none()
        .fill(PANEL)
        .inner_margin(Margin::symmetric(12.0, 10.0))
        .rounding(Rounding::same(10.0))
        .stroke(Stroke::new(1.0, BORDER))
}

/// Returns the color used for a status chip label.
pub fn chip_color(label: &str) -> Color32 {
    match label {
        "Running" => WARNING,
        "Completed" => SUCCESS,
        "Cancelled" => WARNING,
        "Failed" => DANGER,
        _ => Color32::GRAY,
    }
}

/// Maps a container group value to a stable display color.
pub fn group_color(group: u8) -> Color32 {
    if group == 0 {
        return Color32::from_rgb(40, 47, 62);
    }
    let hue = ((group as f32 * 37.0) % 360.0) / 360.0;
    let hsv = egui::ecolor::Hsva::new(hue, 0.72, 0.82, 1.0);
    Color32::from(hsv)
}

/// Returns readable foreground text color for a given background.
pub fn text_on(bg: Color32) -> Color32 {
    let luma = 0.2126 * bg.r() as f32 + 0.7152 * bg.g() as f32 + 0.0722 * bg.b() as f32;
    if luma > 140.0 {
        Color32::BLACK
    } else {
        Color32::WHITE
    }
}
