use std::fs;
use std::path::Path;

fn main() {
    // Ensure dist/flulens.html is always in sync with prototypes/flulens.html
    let candidate_srcs = [
        Path::new("../prototypes/flulens.html"),
        Path::new("prototypes/flulens.html"),
    ];
    let dst_dir = Path::new("dist");
    let dst = dst_dir.join("flulens.html");

    for src in candidate_srcs {
        if src.exists() {
            let _ = fs::create_dir_all(dst_dir);
            let _ = fs::copy(src, &dst);
            println!("cargo:rerun-if-changed={}", src.display());
            break;
        }
    }

    println!("cargo:rerun-if-changed=src/main.rs");
    println!("cargo:rerun-if-changed=tauri.conf.json");
    println!("cargo:rerun-if-changed=capabilities/default.json");
    tauri_build::build();
}
