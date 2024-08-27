use std::{ffi::OsStr,
          fs::File,
          path::Path,
          io::{BufRead, BufReader}};
use flate2::read::MultiGzDecoder;

pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(_why) => panic!("couldn't open file {}", filename),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}