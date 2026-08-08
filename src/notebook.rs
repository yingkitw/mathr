//! Math notebook format (`.mnb`) — a JSON-based notebook with cells of
//! TeX/math expressions and their evaluated results.
//!
//! File structure (`.mnb`):
//! ```json
//! {
//!   "cells": [
//!     { "id": 0, "input": "sin(pi/4)", "output": "0.7071067811865476" },
//!     { "id": 1, "input": "\\frac{1}{2} + \\frac{3}{4}", "output": "1.25" }
//!   ]
//! }
//! ```

use crate::error::{MathError, Result};
use crate::eval::Context;
use crate::repl;
use std::fs;
use std::path::Path;

/// Type of a notebook cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellType {
    /// A math expression or REPL command — evaluated and rendered with KaTeX.
    Math,
    /// A text/markdown cell — displayed as-is, not evaluated.
    Text,
}

impl CellType {
    /// Returns the string tag used in JSON serialization.
    pub fn as_str(&self) -> &'static str {
        match self {
            CellType::Math => "math",
            CellType::Text => "text",
        }
    }

    /// Parses a cell type from its string tag. Defaults to `Math` for unknown values.
    pub fn from_str(s: &str) -> Self {
        match s.trim() {
            "text" => CellType::Text,
            _ => CellType::Math,
        }
    }
}

impl Default for CellType {
    fn default() -> Self {
        CellType::Math
    }
}

/// A single notebook cell.
#[derive(Debug, Clone)]
pub struct NotebookCell {
    pub id: usize,
    pub input: String,
    pub output: String,
    pub cell_type: CellType,
}

/// A math notebook containing multiple cells.
#[derive(Debug, Clone)]
pub struct Notebook {
    pub cells: Vec<NotebookCell>,
}

impl Notebook {
    /// Create an empty notebook.
    pub fn new() -> Self {
        Self { cells: Vec::new() }
    }

    /// Add a new cell with the given input, returning its id.
    pub fn add_cell(&mut self, input: &str) -> usize {
        let id = self.cells.len();
        self.cells.push(NotebookCell {
            id,
            input: input.to_string(),
            output: String::new(),
            cell_type: CellType::Math,
        });
        id
    }

    /// Add a new cell with a specific type, returning its id.
    pub fn add_cell_with_type(&mut self, input: &str, cell_type: CellType) -> usize {
        let id = self.cells.len();
        self.cells.push(NotebookCell {
            id,
            input: input.to_string(),
            output: String::new(),
            cell_type,
        });
        id
    }

    /// Move a cell up by one position. Returns `Ok(())` if moved, or an error
    /// if the cell is already at the top or doesn't exist.
    pub fn move_cell_up(&mut self, id: usize) -> Result<()> {
        if id == 0 {
            return Err(MathError::InvalidArgument("cell already at top".into()));
        }
        if id >= self.cells.len() {
            return Err(MathError::InvalidArgument(format!("cell {} not found", id)));
        }
        self.cells.swap(id, id - 1);
        self.reindex();
        Ok(())
    }

    /// Move a cell down by one position. Returns `Ok(())` if moved, or an error
    /// if the cell is already at the bottom or doesn't exist.
    pub fn move_cell_down(&mut self, id: usize) -> Result<()> {
        if id + 1 >= self.cells.len() {
            return Err(MathError::InvalidArgument("cell already at bottom".into()));
        }
        self.cells.swap(id, id + 1);
        self.reindex();
        Ok(())
    }

    /// Duplicate a cell, inserting the copy immediately after the original.
    /// Returns the id of the new cell.
    pub fn duplicate_cell(&mut self, id: usize) -> Result<usize> {
        if id >= self.cells.len() {
            return Err(MathError::InvalidArgument(format!("cell {} not found", id)));
        }
        let clone = self.cells[id].clone();
        self.cells.insert(id + 1, clone);
        self.reindex();
        Ok(id + 1)
    }

    /// Set the type of a cell.
    pub fn set_cell_type(&mut self, id: usize, cell_type: CellType) -> Result<()> {
        let cell = self
            .cells
            .get_mut(id)
            .ok_or_else(|| MathError::InvalidArgument(format!("cell {} not found", id)))?;
        cell.cell_type = cell_type;
        if cell_type == CellType::Math {
            cell.output.clear();
        }
        Ok(())
    }

    fn reindex(&mut self) {
        for (i, cell) in self.cells.iter_mut().enumerate() {
            cell.id = i;
        }
    }

    /// Evaluate a single cell by id, storing the result. Text cells are skipped.
    /// The context is mutated in-place so that `let` bindings and `fn` definitions
    /// persist across cell evaluations.
    pub fn eval_cell(&mut self, id: usize, ctx: &mut Context) -> Result<()> {
        let cell = self
            .cells
            .get(id)
            .ok_or_else(|| MathError::InvalidArgument(format!("cell {} not found", id)))?;
        if cell.cell_type == CellType::Text {
            // Text cells: the "output" is just the input itself (rendered as text).
            if let Some(cell) = self.cells.get_mut(id) {
                cell.output = cell.input.clone();
            }
            return Ok(());
        }
        let input = cell.input.clone();
        let result = repl::dispatch_with_ctx(&input, ctx)?;
        let output = match result {
            Some(s) => s,
            None => String::new(),
        };
        if let Some(cell) = self.cells.get_mut(id) {
            cell.output = output;
        }
        Ok(())
    }

    /// Evaluate all cells in order, accumulating context state across cells.
    /// Variables and functions defined in earlier cells are available to later cells.
    pub fn eval_all(&mut self, ctx: &mut Context) -> Result<()> {
        let n = self.cells.len();
        for id in 0..n {
            self.eval_cell(id, ctx)?;
        }
        Ok(())
    }

    /// Update a cell's input.
    pub fn set_input(&mut self, id: usize, input: &str) -> Result<()> {
        let cell = self
            .cells
            .get_mut(id)
            .ok_or_else(|| MathError::InvalidArgument(format!("cell {} not found", id)))?;
        cell.input = input.to_string();
        cell.output.clear();
        Ok(())
    }

    /// Remove a cell by id (re-indexes remaining cells).
    pub fn remove_cell(&mut self, id: usize) -> Result<()> {
        if id >= self.cells.len() {
            return Err(MathError::InvalidArgument(format!("cell {} not found", id)));
        }
        self.cells.remove(id);
        self.reindex();
        Ok(())
    }

    /// Load a notebook from a `.mnb` file.
    pub fn load(path: &Path) -> Result<Self> {
        let content = fs::read_to_string(path).map_err(|e| {
            MathError::InvalidArgument(format!("cannot read notebook file: {}", e))
        })?;
        parse_notebook_json(&content)
    }

    /// Save the notebook to a `.mnb` file.
    pub fn save(&self, path: &Path) -> Result<()> {
        let json = self.to_json();
        fs::write(path, json).map_err(|e| {
            MathError::InvalidArgument(format!("cannot write notebook file: {}", e))
        })?;
        Ok(())
    }

    /// Serialize to JSON string.
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\n  \"cells\": [\n");
        for (i, cell) in self.cells.iter().enumerate() {
            s.push_str("    {\n");
            s.push_str(&format!("      \"id\": {},\n", cell.id));
            s.push_str(&format!(
                "      \"input\": {},\n",
                json_escape(&cell.input)
            ));
            s.push_str(&format!(
                "      \"output\": {},\n",
                json_escape(&cell.output)
            ));
            s.push_str(&format!(
                "      \"cell_type\": {}\n",
                json_escape(cell.cell_type.as_str())
            ));
            if i + 1 < self.cells.len() {
                s.push_str("    },\n");
            } else {
                s.push_str("    }\n");
            }
        }
        s.push_str("  ]\n}\n");
        s
    }
}

impl Default for Notebook {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a notebook from a JSON string.
pub fn parse_notebook_json(s: &str) -> Result<Notebook> {
    let mut nb = Notebook::new();
    // Minimal JSON parser for our flat structure
    // Find "cells" array
    let cells_start = s.find("\"cells\"").ok_or_else(|| {
        MathError::InvalidArgument("notebook JSON: missing 'cells' field".into())
    })?;
    let arr_start = s[cells_start..]
        .find('[')
        .ok_or_else(|| MathError::InvalidArgument("notebook JSON: expected '[' after cells".into()))?;
    let arr_start = cells_start + arr_start;
    let arr_end = find_matching_bracket(s, arr_start, '[', ']')?;
    let arr_content = &s[arr_start + 1..arr_end];

    // Split by top-level object boundaries
    let mut depth = 0i32;
    let mut obj_start = None;
    for (i, ch) in arr_content.char_indices() {
        match ch {
            '{' => {
                if depth == 0 {
                    obj_start = Some(i);
                }
                depth += 1;
            }
            '}' => {
                depth -= 1;
                if depth == 0 {
                    if let Some(start) = obj_start {
                        let obj = &arr_content[start..=i];
                        let cell = parse_cell(obj)?;
                        nb.cells.push(cell);
                    }
                    obj_start = None;
                }
            }
            _ => {}
        }
    }

    Ok(nb)
}

fn parse_cell(obj: &str) -> Result<NotebookCell> {
    let id = extract_json_int(obj, "id").unwrap_or(0);
    let input = extract_json_string(obj, "input").unwrap_or_default();
    let output = extract_json_string(obj, "output").unwrap_or_default();
    let cell_type = extract_json_string(obj, "cell_type")
        .map(|s| CellType::from_str(&s))
        .unwrap_or(CellType::Math);
    Ok(NotebookCell { id, input, output, cell_type })
}

fn extract_json_int(obj: &str, key: &str) -> Option<usize> {
    let pattern = format!("\"{}\"", key);
    let pos = obj.find(&pattern)?;
    let rest = &obj[pos + pattern.len()..];
    let colon = rest.find(':')?;
    let rest = &rest[colon + 1..];
    let rest = rest.trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

fn extract_json_string(obj: &str, key: &str) -> Option<String> {
    let pattern = format!("\"{}\"", key);
    let pos = obj.find(&pattern)?;
    let rest = &obj[pos + pattern.len()..];
    let colon = rest.find(':')?;
    let rest = &rest[colon + 1..];
    let rest = rest.trim_start();
    let quote = rest.find('"')?;
    let rest = &rest[quote + 1..];
    // Find unescaped closing quote
    let mut chars = rest.chars().peekable();
    let mut result = String::new();
    let mut escaped = false;
    for ch in chars.by_ref() {
        if escaped {
            match ch {
                'n' => result.push('\n'),
                't' => result.push('\t'),
                'r' => result.push('\r'),
                '"' => result.push('"'),
                '\\' => result.push('\\'),
                '/' => result.push('/'),
                c => result.push(c),
            }
            escaped = false;
        } else if ch == '\\' {
            escaped = true;
        } else if ch == '"' {
            break;
        } else {
            result.push(ch);
        }
    }
    Some(result)
}

fn find_matching_bracket(s: &str, start: usize, open: char, close: char) -> Result<usize> {
    let mut depth = 0i32;
    for (i, ch) in s[start..].char_indices() {
        if ch == open {
            depth += 1;
        } else if ch == close {
            depth -= 1;
            if depth == 0 {
                return Ok(start + i);
            }
        }
    }
    Err(MathError::InvalidArgument(
        "notebook JSON: unmatched bracket".into(),
    ))
}

fn json_escape(s: &str) -> String {
    let mut result = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\t' => result.push_str("\\t"),
            '\r' => result.push_str("\\r"),
            c => result.push(c),
        }
    }
    result.push('"');
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn notebook_create_and_eval() {
        let mut nb = Notebook::new();
        let id = nb.add_cell("sin(pi/4)");
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        let cell = &nb.cells[id];
        assert!(cell.output.contains("0.707"), "output: {}", cell.output);
    }

    #[test]
    fn notebook_eval_tex_input() {
        let mut nb = Notebook::new();
        let id = nb.add_cell(r"\frac{1}{2} + \frac{3}{4}");
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        let cell = &nb.cells[id];
        assert!(cell.output.contains("1.25"), "output: {}", cell.output);
    }

    #[test]
    fn notebook_eval_all() {
        let mut nb = Notebook::new();
        nb.add_cell("1 + 2");
        nb.add_cell("3 * 4");
        nb.add_cell("sin(0)");
        let mut ctx = Context::standard();
        nb.eval_all(&mut ctx).unwrap();
        assert!(nb.cells[0].output.contains("3"));
        assert!(nb.cells[1].output.contains("12"));
        assert!(nb.cells[2].output.contains("0"));
    }

    #[test]
    fn notebook_set_input() {
        let mut nb = Notebook::new();
        let id = nb.add_cell("1 + 1");
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        nb.set_input(id, "2 + 2").unwrap();
        assert_eq!(nb.cells[id].input, "2 + 2");
        assert_eq!(nb.cells[id].output, "");
    }

    #[test]
    fn notebook_remove_cell() {
        let mut nb = Notebook::new();
        nb.add_cell("a");
        nb.add_cell("b");
        nb.add_cell("c");
        nb.remove_cell(1).unwrap();
        assert_eq!(nb.cells.len(), 2);
        assert_eq!(nb.cells[0].id, 0);
        assert_eq!(nb.cells[1].id, 1);
        assert_eq!(nb.cells[1].input, "c");
    }

    #[test]
    fn notebook_json_roundtrip() {
        let mut nb = Notebook::new();
        nb.add_cell("sin(pi/4)");
        nb.add_cell(r"\frac{1}{2}");
        nb.cells[0].output = "0.707...".to_string();
        nb.cells[1].output = "0.5".to_string();

        let json = nb.to_json();
        let nb2 = parse_notebook_json(&json).unwrap();
        assert_eq!(nb2.cells.len(), 2);
        assert_eq!(nb2.cells[0].input, "sin(pi/4)");
        assert_eq!(nb2.cells[0].output, "0.707...");
        assert_eq!(nb2.cells[1].input, r"\frac{1}{2}");
        assert_eq!(nb2.cells[1].output, "0.5");
    }

    #[test]
    fn notebook_json_escape_special() {
        let mut nb = Notebook::new();
        nb.add_cell("a\nb\tc");
        let json = nb.to_json();
        assert!(json.contains("\\n"));
        assert!(json.contains("\\t"));
        let nb2 = parse_notebook_json(&json).unwrap();
        assert_eq!(nb2.cells[0].input, "a\nb\tc");
    }

    #[test]
    fn notebook_save_load_file() {
        let path = std::env::temp_dir().join("mathr_test_notebook.mnb");
        let mut nb = Notebook::new();
        nb.add_cell("1 + 2");
        nb.add_cell("sin(pi/4)");
        nb.cells[0].output = "3".to_string();
        nb.save(&path).unwrap();
        let nb2 = Notebook::load(&path).unwrap();
        assert_eq!(nb2.cells.len(), 2);
        assert_eq!(nb2.cells[0].input, "1 + 2");
        assert_eq!(nb2.cells[0].output, "3");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn notebook_eval_diff_command() {
        let mut nb = Notebook::new();
        let id = nb.add_cell("diff x^3");
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        let cell = &nb.cells[id];
        assert!(
            cell.output.contains("3") && cell.output.contains("x"),
            "output should contain derivative: {}",
            cell.output
        );
    }

    #[test]
    fn notebook_eval_solve_command() {
        let mut nb = Notebook::new();
        let id = nb.add_cell("solve x^2 - 4");
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        let cell = &nb.cells[id];
        assert!(
            cell.output.contains("2") || cell.output.contains("root"),
            "output should contain root: {}",
            cell.output
        );
    }

    #[test]
    fn notebook_parse_empty_cells() {
        let json = r#"{"cells": []}"#;
        let nb = parse_notebook_json(json).unwrap();
        assert_eq!(nb.cells.len(), 0);
    }

    #[test]
    fn notebook_parse_bad_json() {
        assert!(parse_notebook_json(r#"{"foo": "bar"}"#).is_err());
        assert!(parse_notebook_json(r#"{"cells": }"#).is_err());
    }

    #[test]
    fn notebook_move_cell_up() {
        let mut nb = Notebook::new();
        nb.add_cell("a");
        nb.add_cell("b");
        nb.add_cell("c");
        nb.move_cell_up(2).unwrap(); // c moves to position 1
        assert_eq!(nb.cells[0].input, "a");
        assert_eq!(nb.cells[1].input, "c");
        assert_eq!(nb.cells[2].input, "b");
        assert_eq!(nb.cells[0].id, 0);
        assert_eq!(nb.cells[1].id, 1);
        assert_eq!(nb.cells[2].id, 2);
    }

    #[test]
    fn notebook_move_cell_down() {
        let mut nb = Notebook::new();
        nb.add_cell("a");
        nb.add_cell("b");
        nb.add_cell("c");
        nb.move_cell_down(0).unwrap(); // a moves to position 1
        assert_eq!(nb.cells[0].input, "b");
        assert_eq!(nb.cells[1].input, "a");
        assert_eq!(nb.cells[2].input, "c");
    }

    #[test]
    fn notebook_move_cell_boundaries() {
        let mut nb = Notebook::new();
        nb.add_cell("a");
        nb.add_cell("b");
        assert!(nb.move_cell_up(0).is_err()); // already at top
        assert!(nb.move_cell_down(1).is_err()); // already at bottom
        assert!(nb.move_cell_up(99).is_err()); // doesn't exist
    }

    #[test]
    fn notebook_duplicate_cell() {
        let mut nb = Notebook::new();
        nb.add_cell("a");
        nb.add_cell("b");
        let new_id = nb.duplicate_cell(0).unwrap();
        assert_eq!(nb.cells.len(), 3);
        assert_eq!(nb.cells[0].input, "a");
        assert_eq!(nb.cells[1].input, "a"); // duplicate
        assert_eq!(nb.cells[2].input, "b");
        assert_eq!(new_id, 1);
    }

    #[test]
    fn notebook_cell_types() {
        let mut nb = Notebook::new();
        let math_id = nb.add_cell("sin(pi/4)");
        let text_id = nb.add_cell_with_type("# My Heading", CellType::Text);
        assert_eq!(nb.cells[math_id].cell_type, CellType::Math);
        assert_eq!(nb.cells[text_id].cell_type, CellType::Text);
    }

    #[test]
    fn notebook_text_cell_eval() {
        let mut nb = Notebook::new();
        let id = nb.add_cell_with_type("This is documentation", CellType::Text);
        nb.eval_cell(id, &mut Context::standard()).unwrap();
        // Text cells: output is the input itself
        assert_eq!(nb.cells[id].output, "This is documentation");
    }

    #[test]
    fn notebook_set_cell_type() {
        let mut nb = Notebook::new();
        let id = nb.add_cell("1 + 2");
        nb.cells[id].output = "3".to_string();
        nb.set_cell_type(id, CellType::Text).unwrap();
        assert_eq!(nb.cells[id].cell_type, CellType::Text);
        // Switching back to Math clears output
        nb.set_cell_type(id, CellType::Math).unwrap();
        assert_eq!(nb.cells[id].output, "");
    }

    #[test]
    fn notebook_json_with_cell_types() {
        let mut nb = Notebook::new();
        nb.add_cell("sin(pi/4)");
        nb.add_cell_with_type("# Notes", CellType::Text);
        let json = nb.to_json();
        assert!(json.contains("\"cell_type\""));
        assert!(json.contains("math"));
        assert!(json.contains("text"));
        let nb2 = parse_notebook_json(&json).unwrap();
        assert_eq!(nb2.cells[0].cell_type, CellType::Math);
        assert_eq!(nb2.cells[1].cell_type, CellType::Text);
    }

    #[test]
    fn notebook_json_backward_compatible_no_cell_type() {
        // Old .mnb files without cell_type should default to Math
        let json = r#"{"cells": [{"id": 0, "input": "1+2", "output": "3"}]}"#;
        let nb = parse_notebook_json(json).unwrap();
        assert_eq!(nb.cells[0].cell_type, CellType::Math);
    }

    #[test]
    fn notebook_eval_all_shared_context() {
        let mut nb = Notebook::new();
        nb.add_cell("let x = 5");
        nb.add_cell("x * 2");
        let mut ctx = Context::standard();
        nb.eval_all(&mut ctx).unwrap();
        // Second cell should see x = 5 from the first cell
        assert!(nb.cells[1].output.contains("10"), "output: {}", nb.cells[1].output);
    }
}
