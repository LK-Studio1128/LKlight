use lklight::simulator::{parse_method, read_setup_from_file, simulate};
use std::env;
use std::path::Path;
use std::thread;

const STACK_SIZE: usize = 8 * 1024 * 1024;

fn main() {
    let child = thread::Builder::new().stack_size(STACK_SIZE).spawn(run).unwrap();
    child.join().unwrap();
}

fn run() {
    env_logger::init();
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} setup_filename swarm_filename steps method", args[0]);
        return;
    }
    let steps: u32 = match args[3].parse() {
        Ok(n) => n,
        Err(_) => {
            eprintln!("Error: steps must be a number");
            return;
        }
    };
    let method = match parse_method(&args[4]) {
        Some(m) => m,
        None => {
            eprintln!("Unsupported method. Supported: ddna dfire fastdfire dfire2 dna mj3h pisa pydock cpydock sipper tobi vdw");
            return;
        }
    };
    let setup = match read_setup_from_file(&args[1]) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error reading setup file: {:?}", e);
            return;
        }
    };
    let sim_path = Path::new(&args[1]).parent().unwrap();
    simulate(sim_path.to_str().unwrap(), &setup, &args[2], steps, method);
}
