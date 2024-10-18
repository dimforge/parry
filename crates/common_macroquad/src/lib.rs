#[cfg(feature = "dim2")]
pub mod dim2;
#[cfg(feature = "dim3")]
pub mod dim3;

use macroquad::color::WHITE;

/// Draws a text in the top left corner of the screen.
///
/// This uses a hardcoded position, size, color.
pub fn easy_draw_text(text: &str) {
    macroquad::text::draw_text(text, 10.0, 48.0 + 18.0, 30.0, WHITE);
}
