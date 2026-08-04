// FluLens desktop shell.
//
// The whole job of this file is to hand the page the same three operations the
// browser's File object gives it — read the text, read a byte range, ask the
// size — over a directory the user picked, and to remember which directory that
// was. Everything else stays in flulens.html.
//
// The byte-range read is not a convenience. A BAM is hundreds of megabytes and
// the app queries a few kilobytes of it per codon; reading whole files would
// make the pileup unusable. In the browser that range comes from File.slice, in
// dev from an HTTP Range header, and here from a seek.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use tauri::Manager;

/// Directories that are pipeline scratch, never results.
///
/// Nextflow's work/ sits inside the run directory and holds a staged copy of
/// nearly everything — on the swine run, 872 copies of reference.fa and 429
/// VCFs. Every lookup in the app takes the first path that matches, so those
/// copies do not sit inert, they compete: reference_gtf was read twice and PB2
/// came out as 1518 codons instead of 759. The browser path filters these in
/// indexFiles(); this is the same rule, applied before the paths are ever sent.
fn is_scratch(rel: &str) -> bool {
    rel.split('/').any(|p| p == "work" || p == ".nextflow")
}

fn walk(root: &Path, dir: &Path, out: &mut Vec<String>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return, // unreadable subdirectory is not fatal — skip it
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let rel = match path.strip_prefix(root) {
            Ok(r) => r.to_string_lossy().replace('\\', "/"),
            Err(_) => continue,
        };
        if is_scratch(&rel) {
            continue;
        }
        match entry.file_type() {
            Ok(t) if t.is_dir() => walk(root, &path, out),
            Ok(t) if t.is_file() => out.push(rel),
            _ => {}
        }
    }
}

/// Every readable file under `root`, as run-relative paths, scratch excluded.
#[tauri::command]
fn list_files(root: String) -> Result<Vec<String>, String> {
    let rootp = PathBuf::from(&root);
    if !rootp.is_dir() {
        return Err(format!("not a directory: {root}"));
    }
    let mut out = Vec::new();
    walk(&rootp, &rootp, &mut out);
    Ok(out)
}

#[tauri::command]
fn read_text(root: String, rel: String) -> Result<String, String> {
    let p = safe_join(&root, &rel)?;
    fs::read(&p)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .map_err(|e| format!("{}: {e}", p.display()))
}

/// Bytes [start, end) — the seek that makes BAM queries cheap.
#[tauri::command]
fn read_bytes(root: String, rel: String, start: u64, end: u64) -> Result<Vec<u8>, String> {
    let p = safe_join(&root, &rel)?;
    if end <= start {
        return Ok(Vec::new());
    }
    let mut f = fs::File::open(&p).map_err(|e| format!("{}: {e}", p.display()))?;
    f.seek(SeekFrom::Start(start)).map_err(|e| e.to_string())?;
    let mut buf = vec![0u8; (end - start) as usize];
    // read_exact would fail at EOF; a short final chunk is normal and fine.
    let mut filled = 0usize;
    loop {
        match f.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => {
                filled += n;
                if filled >= buf.len() {
                    break;
                }
            }
            Err(e) => return Err(e.to_string()),
        }
    }
    buf.truncate(filled);
    Ok(buf)
}

#[tauri::command]
fn file_size(root: String, rel: String) -> Result<u64, String> {
    let p = safe_join(&root, &rel)?;
    fs::metadata(&p).map(|m| m.len()).map_err(|e| e.to_string())
}

/// Join and verify the result is still inside `root`.
///
/// The page supplies `rel`, and the page renders file paths that came from a
/// results directory. Nothing there should be able to address `../../` its way
/// out of the run, whether by malice or by a malformed path in a table.
fn safe_join(root: &str, rel: &str) -> Result<PathBuf, String> {
    let rootp = fs::canonicalize(root).map_err(|e| e.to_string())?;
    let joined = rootp.join(rel);
    let real = fs::canonicalize(&joined).map_err(|e| format!("{}: {e}", joined.display()))?;
    if !real.starts_with(&rootp) {
        return Err(format!("path escapes the run directory: {rel}"));
    }
    Ok(real)
}

/// The last opened run, so the app does not ask again every launch. This is the
/// single thing the browser build cannot do: a File handle dies with the tab.
#[tauri::command]
fn last_run(app: tauri::AppHandle) -> Option<String> {
    let p = app.path().app_config_dir().ok()?.join("last_run.txt");
    fs::read_to_string(p).ok().map(|s| s.trim().to_string()).filter(|s| !s.is_empty())
}

#[tauri::command]
fn set_last_run(app: tauri::AppHandle, path: String) -> Result<(), String> {
    let dir = app.path().app_config_dir().map_err(|e| e.to_string())?;
    fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    fs::write(dir.join("last_run.txt"), path).map_err(|e| e.to_string())
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            list_files, read_text, read_bytes, file_size, last_run, set_last_run
        ])
        .run(tauri::generate_context!())
        .expect("error while running FluLens");
}
