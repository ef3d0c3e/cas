//! Loading libgiac at run time and calling its C entry points: the only
//! `unsafe` code of the crate, and of the workspace, allowed here because
//! there is no other way to call a C library loaded at run time. giac's
//! headers declare, in `gen.h`:
//!
//! ```c
//! extern "C" const char * caseval(const char *);
//! extern "C" const char * nws_caseval(const char * s);
//! ```
//!
//! Each returns a static buffer that the next call overwrites, and giac
//! keeps global state, so every call takes [`LOCK`] and copies the answer
//! before releasing it.

use std::ffi::{CStr, CString, c_char, c_void};
use std::sync::Mutex;

/// The signature of `caseval` and `nws_caseval`.
type Entry = unsafe extern "C" fn(*const c_char) -> *const c_char;

/// Held across each call and the copy of its answer.
static LOCK: Mutex<()> = Mutex::new(());

/// The two entry points of a loaded libgiac. The library stays loaded for
/// the rest of the process: its handle is never closed, so the entry points
/// stay valid.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Library {
    caseval: Entry,
    nws_caseval: Entry,
}

/// The text of `dlerror`, or `what` when it has none.
fn last_error(what: &str) -> String {
    // SAFETY: `dlerror` returns null or a NUL-terminated string valid until
    // the next `dl*` call on this thread, which comes after the copy.
    let message = unsafe { libc::dlerror() };
    if message.is_null() {
        return what.to_owned();
    }
    // SAFETY: as above, non-null.
    unsafe { CStr::from_ptr(message) }
        .to_string_lossy()
        .into_owned()
}

impl Library {
    /// The library at `path` (a file name searches the usual paths).
    pub(crate) fn load(path: &str) -> Result<Library, String> {
        let name = CString::new(path).map_err(|_| format!("{path:?}: a NUL in the path"))?;
        let _held = LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        // SAFETY: `name` is a valid NUL-terminated string. Loading runs the
        // library's initialisers, giac's, which is what using it takes.
        let handle = unsafe { libc::dlopen(name.as_ptr(), libc::RTLD_NOW | libc::RTLD_LOCAL) };
        if handle.is_null() {
            return Err(last_error("dlopen failed"));
        }
        let symbol = |symbol: &CStr| -> Result<Entry, String> {
            // SAFETY: `handle` is a live handle that is never closed, and
            // `symbol` a valid NUL-terminated string.
            let address: *mut c_void = unsafe { libc::dlsym(handle, symbol.as_ptr()) };
            if address.is_null() {
                return Err(last_error("dlsym failed"));
            }
            // SAFETY: giac's `gen.h` declares both symbols as functions of
            // this C signature, and the address is that of the function.
            Ok(unsafe { std::mem::transmute::<*mut c_void, Entry>(address) })
        };
        Ok(Library {
            caseval: symbol(c"caseval")?,
            nws_caseval: symbol(c"nws_caseval")?,
        })
    }

    /// `caseval(input)`, or `nws_caseval(input)` when `fast`; `None` for an
    /// input with a NUL, or a null answer.
    pub(crate) fn call(&self, input: &str, fast: bool) -> Option<String> {
        let input = CString::new(input).ok()?;
        let entry = if fast { self.nws_caseval } else { self.caseval };
        let _held = LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        // SAFETY: `entry` is one of giac's two entry points, of the declared
        // signature, and the library is never unloaded; `input` is a valid
        // NUL-terminated string that outlives the call. The lock keeps other
        // threads out of giac's global state.
        let answer = unsafe { entry(input.as_ptr()) };
        if answer.is_null() {
            return None;
        }
        // SAFETY: a non-null answer is a NUL-terminated static buffer, valid
        // until the next call, which the lock holds off until the copy.
        Some(
            unsafe { CStr::from_ptr(answer) }
                .to_string_lossy()
                .into_owned(),
        )
    }
}
