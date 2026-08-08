//! Plotting via the [`plotters`] crate. Output is rendered to PNG files.
//!
//! Three plot kinds are supported:
//! - [`plot_function`]: 2D line plot of an `f(x)` over a range
//! - [`plot_scatter`]: scatter or stem plot, useful for FFT magnitude / power
//! - [`plot_multi`]: several functions overlaid on a single chart

use crate::error::{MathError, Result};
use crate::expr::Expr;
use plotters::prelude::*;
use std::path::Path;

/// Render a plot to a PNG file, read the bytes into memory, and delete the temp file.
/// Used by the notebook server for inline image display.
fn plot_to_bytes<F>(draw: F) -> Result<Vec<u8>>
where
    F: FnOnce(&Path) -> Result<()>,
{
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let id = COUNTER.fetch_add(1, Ordering::SeqCst);
    let path = std::env::temp_dir().join(format!("mathr_plot_{}_{}.png", std::process::id(), id));
    draw(&path)?;
    let bytes = std::fs::read(&path)
        .map_err(|e| MathError::Plot(format!("cannot read temp plot: {}", e)))?;
    let _ = std::fs::remove_file(&path);
    Ok(bytes)
}

/// Plot a single expression to PNG bytes in memory (for inline notebook display).
pub fn plot_function_to_bytes(
    expr: &Expr,
    x_var: &str,
    x_min: f64,
    x_max: f64,
    samples: usize,
    title: &str,
) -> Result<Vec<u8>> {
    plot_to_bytes(|path| {
        plot_function(path, expr, x_var, x_min, x_max, samples, title)
    })
}

/// Plot multiple series to PNG bytes in memory.
pub fn plot_multi_to_bytes(
    series: &[(String, Expr, &str)],
    x_min: f64,
    x_max: f64,
    samples: usize,
    title: &str,
) -> Result<Vec<u8>> {
    plot_to_bytes(|path| {
        plot_multi(path, series, x_min, x_max, samples, title)
    })
}

/// Plot a scatter to PNG bytes in memory.
pub fn plot_scatter_to_bytes(
    points: &[(f64, f64)],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<Vec<u8>> {
    plot_to_bytes(|path| {
        plot_scatter(path, points, title, x_label, y_label)
    })
}

/// Convert any plotters error to a MathError so the `?` operator can
/// propagate drawing-area failures back through our `Result` type.
fn p_err<E: std::fmt::Display>(e: E) -> MathError {
    MathError::Plot(e.to_string())
}

/// Plot a single expression `f(x_var)` over `[x_min, x_max]` to a PNG file.
/// `samples` is the number of points to evaluate.
pub fn plot_function<P: AsRef<Path>>(
    path: P,
    expr: &Expr,
    x_var: &str,
    x_min: f64,
    x_max: f64,
    samples: usize,
    title: &str,
) -> Result<()> {
    if samples < 2 {
        return Err(MathError::InvalidArgument("plot needs at least 2 samples".into()));
    }
    if !(x_min.is_finite() && x_max.is_finite() && x_min < x_max) {
        return Err(MathError::InvalidArgument(format!(
            "plot: bad x range [{}, {}]",
            x_min, x_max
        )));
    }

    let ctx = crate::eval::Context::standard();
    let mut points: Vec<(f64, f64)> = Vec::with_capacity(samples);
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for i in 0..samples {
        let x = x_min + (x_max - x_min) * i as f64 / (samples - 1) as f64;
        let mut cx = ctx.clone();
        cx.set(x_var, x);
        match crate::eval::eval(expr, &cx) {
            Ok(y) if y.is_finite() => {
                points.push((x, y));
                if y < y_min {
                    y_min = y;
                }
                if y > y_max {
                    y_max = y;
                }
            }
            _ => {} // skip asymptotes / domain errors
        }
    }

    if points.is_empty() {
        return Err(MathError::Plot(
            "no valid samples to plot (function undefined over entire range?)".into(),
        ));
    }

    let pad = (y_max - y_min).max(1e-6) * 0.1;
    y_min -= pad;
    y_max += pad;

    let root = BitMapBackend::new(&path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE).map_err(p_err)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(p_err)?;
    chart
        .configure_mesh()
        .x_label_formatter(&|x| format!("{:.2}", x))
        .y_label_formatter(&|y| format!("{:.2}", y))
        .draw()
        .map_err(p_err)?;

    chart
        .draw_series(LineSeries::new(points, &RED))
        .map_err(p_err)?
        .label(title)
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .map_err(p_err)?;
    Ok(())
}

/// Plot several expressions on the same axes. Useful for comparing functions
/// or overlaying a derivative on its source.
pub fn plot_multi<P: AsRef<Path>>(
    path: P,
    series: &[(String, Expr, &str)],
    x_min: f64,
    x_max: f64,
    samples: usize,
    title: &str,
) -> Result<()> {
    if !(x_min.is_finite() && x_max.is_finite() && x_min < x_max) {
        return Err(MathError::InvalidArgument(format!(
            "plot_multi: bad x range [{}, {}]",
            x_min, x_max
        )));
    }
    if series.is_empty() {
        return Err(MathError::InvalidArgument("plot_multi needs at least one series".into()));
    }
    let ctx = crate::eval::Context::standard();

    // Sample every series.
    let mut data: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for (_, expr, var) in series {
        let mut pts = Vec::with_capacity(samples);
        for i in 0..samples {
            let x = x_min + (x_max - x_min) * i as f64 / (samples - 1) as f64;
            let mut cx = ctx.clone();
            cx.set(*var, x);
            if let Ok(y) = crate::eval::eval(expr, &cx) {
                if y.is_finite() {
                    pts.push((x, y));
                    y_min = y_min.min(y);
                    y_max = y_max.max(y);
                }
            }
        }
        data.push(pts);
    }
    if !y_min.is_finite() || !y_max.is_finite() {
        return Err(MathError::Plot("no valid samples for any series".into()));
    }
    let pad = (y_max - y_min).max(1e-6) * 0.1;
    y_min -= pad;
    y_max += pad;

    let root = BitMapBackend::new(&path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE).map_err(p_err)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(p_err)?;
    chart.configure_mesh().draw().map_err(p_err)?;

    // Pick a fixed palette in order, indexable by position.
    let palette: [RGBColor; 6] = [RED, BLUE, GREEN, MAGENTA, CYAN, BLACK];
    for (idx, item) in series.iter().enumerate() {
        let name = &item.0;
        let points = &data[idx];
        let color = palette[idx % palette.len()];
        let legend_color = color.clone();
        chart
            .draw_series(LineSeries::new(points.iter().copied(), color.clone()))
            .map_err(p_err)?
            .label(name.as_str())
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], legend_color.clone())
            });
    }
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .map_err(p_err)?;
    Ok(())
}

/// Plot a scatter series — typically used for FFT magnitudes.
pub fn plot_scatter<P: AsRef<Path>>(
    path: P,
    points: &[(f64, f64)],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    if points.is_empty() {
        return Err(MathError::InvalidArgument(
            "plot_scatter needs at least one point".into(),
        ));
    }
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for &(x, y) in points {
        x_min = x_min.min(x);
        x_max = x_max.max(x);
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }
    let x_pad = (x_max - x_min).max(1e-6) * 0.05;
    let y_pad = (y_max - y_min).max(1e-6) * 0.10;
    x_min -= x_pad;
    x_max += x_pad;
    y_min -= y_pad;
    y_max += y_pad;

    let root = BitMapBackend::new(&path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE).map_err(p_err)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(p_err)?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()
        .map_err(p_err)?;
    chart
        .draw_series(
            points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 3, BLUE.filled())),
        )
        .map_err(p_err)?
        .label(title)
        .legend(|(x, y)| Circle::new((x + 10, y), 3, BLUE.filled()));
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .map_err(p_err)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    #[test]
    fn plot_function_to_bytes_produces_png() {
        let expr = Parser::parse("sin(x)").unwrap();
        let bytes = plot_function_to_bytes(&expr, "x", 0.0, 3.14159, 100, "y = sin(x)").unwrap();
        assert!(bytes.len() > 100, "PNG should have content");
        assert_eq!(&bytes[..8], &[137, 80, 78, 71, 13, 10, 26, 10], "should be PNG");
    }

    #[test]
    fn plot_function_to_bytes_bad_range() {
        let expr = Parser::parse("x").unwrap();
        let result = plot_function_to_bytes(&expr, "x", 5.0, 1.0, 100, "y = x");
        assert!(result.is_err());
    }

    #[test]
    fn plot_scatter_to_bytes_produces_png() {
        let points = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];
        let bytes = plot_scatter_to_bytes(&points, "quadratic", "x", "y").unwrap();
        assert!(bytes.len() > 100);
        assert_eq!(&bytes[..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    }
}