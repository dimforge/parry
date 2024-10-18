mod plane_intersection;
use fork::{fork, Fork};

fn main() {
    // NOTE: we can't use threads via macroquad: it supports only 1 window ; so I'm forking it.
    // This will be problematic to support wasm...
    // ...but! maybe not entirely difficult to make them target a specific canvas
    // This would still need a custom communication layer between both windows/canvas,
    // which would be different between web and native
    match fork() {
        Ok(Fork::Parent(child)) => {
            plane_intersection::main();
            println!(
                "Continuing execution in parent process, new child has pid: {}",
                child
            );
        }
        Ok(Fork::Child) => {
            plane_intersection::main();
            println!("I'm a new child process")
        }
        Err(_) => println!("Fork failed"),
    }
}
