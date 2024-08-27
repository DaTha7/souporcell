pub fn relative_index_search(normalized_logs: &Vec<f32>,
                             log_value_to_find: f32) -> Option<usize> {

    let length = normalized_logs.len();
    let mut mid = length / 2;
    let mut hind = length - 1;
    let mut lind = 0;
    let mut current = normalized_logs[mid];

    while lind <= hind {

        match current.total_cmp(&log_value_to_find) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => lind = mid + 1,
            std::cmp::Ordering::Greater => hind = mid - 1,
        }

        mid = (hind + lind) / 2;
        current = normalized_logs[mid];
    }

    Some(lind)
    // Some(hind)
}