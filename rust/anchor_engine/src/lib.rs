use std::collections::BTreeSet;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;

fn split_segments(motif: &str, max_mismatches: usize) -> Vec<(usize, &str)> {
    let segment_count = max_mismatches + 1;
    let base_size = motif.len() / segment_count;
    let remainder = motif.len() % segment_count;
    let mut segments = Vec::with_capacity(segment_count);
    let mut start = 0usize;
    for idx in 0..segment_count {
        let seg_len = base_size + usize::from(idx < remainder);
        if seg_len == 0 {
            continue;
        }
        let end = start + seg_len;
        segments.push((start, &motif[start..end]));
        start = end;
    }
    segments
}

fn bounded_hamming(left: &[u8], right: &[u8], max_mismatches: usize) -> Option<usize> {
    let mut mismatches = 0usize;
    for (lch, rch) in left.iter().zip(right.iter()) {
        if lch != rch {
            mismatches += 1;
            if mismatches > max_mismatches {
                return None;
            }
        }
    }
    Some(mismatches)
}

fn find_exact_hits(sequence: &str, motif: &str) -> String {
    let mut out = String::new();
    for (idx, _) in sequence.match_indices(motif) {
        if !out.is_empty() {
            out.push(',');
        }
        out.push_str(&format!("{}:0", idx));
    }
    out
}

fn find_approximate_hits(sequence: &str, motif: &str, max_mismatches: usize) -> String {
    let segments = split_segments(motif, max_mismatches);
    let seq_len = sequence.len();
    let motif_len = motif.len();
    let mut candidate_starts = BTreeSet::new();
    for (segment_offset, segment) in segments {
        for (idx, _) in sequence.match_indices(segment) {
            if idx >= segment_offset {
                let start = idx - segment_offset;
                if start + motif_len <= seq_len {
                    candidate_starts.insert(start);
                }
            }
        }
    }

    let seq_bytes = sequence.as_bytes();
    let motif_bytes = motif.as_bytes();
    let mut out = String::new();
    for start in candidate_starts {
        let end = start + motif_len;
        if let Some(mismatches) = bounded_hamming(&seq_bytes[start..end], motif_bytes, max_mismatches) {
            if !out.is_empty() {
                out.push(',');
            }
            out.push_str(&format!("{}:{}", start, mismatches));
        }
    }
    out
}

#[no_mangle]
pub extern "C" fn anchorscope_find_fixed_hits(
    sequence: *const c_char,
    motif: *const c_char,
    max_mismatches: usize,
) -> *mut c_char {
    if sequence.is_null() || motif.is_null() {
        return std::ptr::null_mut();
    }

    let sequence = unsafe { CStr::from_ptr(sequence) }.to_string_lossy().to_uppercase();
    let motif = unsafe { CStr::from_ptr(motif) }.to_string_lossy().to_uppercase();

    let result = if max_mismatches == 0 {
        find_exact_hits(&sequence, &motif)
    } else {
        find_approximate_hits(&sequence, &motif, max_mismatches)
    };

    match CString::new(result) {
        Ok(c_string) => c_string.into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn anchorscope_free_string(ptr: *mut c_char) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        let _ = CString::from_raw(ptr);
    }
}
