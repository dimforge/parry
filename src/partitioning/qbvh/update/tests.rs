use rand::{rngs::StdRng, Rng, SeedableRng};
use std::borrow::Cow;

use crate::{
    bounding_volume::Aabb,
    math::{Point, Real},
    partitioning::{GenericQbvh, Qbvh},
    utils::DefaultStorage,
};

use super::QbvhUpdateWorkspace;

#[test]
fn test_case_qbvh_1() {
    //  thread 'partitioning::qbvh::update::tests::test_case_qbvh' panicked at 'assertion failed: `(left == right)`
    //  left: `5`,
    //  right: `6`', crates\parry2d\../../src\partitioning\qbvh\update.rs:252:9
    test_qbvh_random_operations(0x23fc68663e15e9e2, 100, false, None, 0.55, 200);
}

#[test]
fn test_case_qbvh_2() {
    // thread 'partitioning::qbvh::update::tests::test_case_qbvh' panicked at 'failed for 977', crates\parry2d\../../src\partitioning\qbvh\update.rs:203:21
    test_qbvh_random_operations(0xd758a3c214dc3866, 2_406, false, None, 0.55, 1000);
}

#[test]
fn test_case_qbvh_3() {
    // thread 'partitioning::qbvh::update::tests::test_case_qbvh_3' panicked at 'attempt to add with overflow', crates\parry2d\../../src\partitioning\qbvh\update.rs:363:54
    test_qbvh_random_operations(0x93bcaea3b92a9bfe, 100, false, None, 0.53, 100);
}

// // uncomment to search for a failed seed
// #[test]
// fn test_case_qbvh() {
//     // do this in multiple thread
//     let found = std::sync::Arc::new(core::sync::atomic::AtomicBool::new(false));
//     for _ in 0..std::thread::available_parallelism().unwrap().into() {
//         let found = found.clone();

//         let _handle = std::thread::spawn(move || loop {
//             let seed = rand::random();
//             let result = std::panic::catch_unwind(|| {
//                 test_qbvh_random_operations(seed, 10000, false, None, 0.53, 1000);
//             });

//             match result {
//                 Ok(()) => {}
//                 Err(err) => {
//                     found.store(true, core::sync::atomic::Ordering::Relaxed);
//                     let message = if let Some(msg) = err.downcast_ref::<&str>() {
//                         Cow::<str>::from(*msg)
//                     } else if let Some(msg) = err.downcast_ref::<String>() {
//                         Cow::<str>::from(msg.clone())
//                     } else {
//                         Cow::<str>::from("Unknown error")
//                     };
//                     panic!("Test failed with seed: 0x{seed:x} because {message}");
//                 }
//             }
//         });
//     }
//     loop {
//         if found.load(core::sync::atomic::Ordering::Relaxed) {
//             panic!("Found a failed seed");
//         }
//     }
// }

fn test_qbvh_random_operations(
    seed: u64,
    num_iteration: usize,
    debug_print: bool,
    debug_specific_iteration: Option<usize>,
    grow_ratio: f32,
    max_aabbs_count: usize,
) {
    let mut rng = StdRng::seed_from_u64(seed);

    let mut qbvh = QbvhTester::new();
    let mut added_aabb_indices = vec![];

    for count in 0..num_iteration {
        let x = rng.gen_range(0.0..1.0);
        if x < grow_ratio && added_aabb_indices.len() < max_aabbs_count {
            let aabb: Aabb = generate_random_aabb(&mut rng);
            if debug_print {
                #[cfg(feature = "dim3")]
                println!(
                    "{:5}, ADD, mins: [{:>4.0}, {:>4.0}, {:>4.0}], maxs: [{:>4.0}, {:>4.0}, {:>4.0}]",
                    count, aabb.mins.x, aabb.mins.y, aabb.mins.z, aabb.maxs.x, aabb.maxs.y, aabb.maxs.z,
                );
                #[cfg(feature = "dim2")]
                println!(
                    "{:5}, ADD, mins: [{:>4.0}, {:>4.0}], maxs: [{:>4.0}, {:>4.0}]",
                    count, aabb.mins.x, aabb.mins.y, aabb.maxs.x, aabb.maxs.y,
                );
            }
            let index;
            if let Some(iteration) = debug_specific_iteration {
                if count == iteration {
                    {
                        println!("before insert");
                        qbvh.print_tree();
                    }
                    {
                        index = qbvh.aabbs.insert(aabb);
                        qbvh.qbvh.pre_update_or_insert(index);
                        println!("before refit");
                        qbvh.print_tree();
                    }
                    {
                        let _count = qbvh
                            .qbvh
                            .refit(0.0, &mut qbvh.workspace, |index| qbvh.aabbs[*index]);
                        println!("before rebalance");
                        qbvh.print_tree();
                    }
                    {
                        qbvh.qbvh.rebalance(0.0, &mut qbvh.workspace);
                        println!("after rebalance");
                        qbvh.print_tree();
                    }
                    {
                        let cloned = qbvh.get_clear_and_rebuild();
                        cloned.check_topology();
                        println!("clean");
                        cloned.print_tree();
                    }
                } else {
                    index = qbvh.add_aabb(aabb);
                }
            } else {
                index = qbvh.add_aabb(aabb);
            };
            added_aabb_indices.push(index);
            qbvh.check_topology();
        } else {
            // remove aabb
            if added_aabb_indices.is_empty() {
                continue;
            }
            let aabb_index =
                added_aabb_indices.swap_remove(rng.gen_range(0..added_aabb_indices.len()));
            if debug_print {
                let aabb = qbvh.aabbs[aabb_index];
                #[cfg(feature = "dim3")]
                println!(
                    "{:5}, RMV, mins: [{:>4.0}, {:>4.0}, {:>4.0}], maxs: [{:>4.0}, {:>4.0}, {:>4.0}], index: {:<5}",
                    count, aabb.mins.x, aabb.mins.y, aabb.mins.z, aabb.maxs.x, aabb.maxs.y, aabb.maxs.z, aabb_index
                );
                #[cfg(feature = "dim2")]
                println!(
                    "{:5}, RMV, mins: [{:>4.0}, {:>4.0}], maxs: [{:>4.0}, {:>4.0}], index: {:<5}",
                    count, aabb.mins.x, aabb.mins.y, aabb.maxs.x, aabb.maxs.y, aabb_index
                );
            }
            if let Some(iteration) = debug_specific_iteration {
                if count == iteration {
                    {
                        println!("before insert");
                        qbvh.print_tree();
                    }
                    {
                        let _aabb = qbvh.aabbs.remove(aabb_index);
                        let _removed = qbvh.qbvh.remove(aabb_index);
                        println!("before refit");
                        qbvh.print_tree();
                    }
                    {
                        let _count = qbvh
                            .qbvh
                            .refit(0.0, &mut qbvh.workspace, |index| qbvh.aabbs[*index]);
                        println!("before rebalance");
                        qbvh.print_tree();
                    }
                    {
                        qbvh.qbvh.rebalance(0.0, &mut qbvh.workspace);
                        println!("after rebalance");
                        qbvh.print_tree();
                    }
                    {
                        let cloned = qbvh.get_clear_and_rebuild();
                        println!("clean");
                        cloned.print_tree();
                    }
                } else {
                    qbvh.remove_aabb(aabb_index);
                }
            } else {
                qbvh.remove_aabb(aabb_index);
            }
            qbvh.check_topology();
        }
    }
}

pub struct QbvhTester {
    pub qbvh: Qbvh<usize>,
    pub workspace: QbvhUpdateWorkspace,
    pub aabbs: slab::Slab<Aabb>,
}

impl QbvhTester {
    fn new() -> Self {
        Self {
            qbvh: Qbvh::new(),
            workspace: QbvhUpdateWorkspace::default(),
            aabbs: slab::Slab::new(),
        }
    }

    fn get_clear_and_rebuild(&self) -> QbvhTester {
        let mut qbvh = Qbvh::new();
        let workspace = QbvhUpdateWorkspace::default();
        let aabbs = self.aabbs.clone();
        qbvh.clear_and_rebuild(aabbs.iter().map(|(index, aabb)| (index, aabb.clone())), 0.0);
        QbvhTester {
            qbvh,
            workspace,
            aabbs,
        }
    }

    fn add_aabb(&mut self, aabb: Aabb) -> usize {
        let index = self.aabbs.insert(aabb);
        self.qbvh.pre_update_or_insert(index);
        let _count = self
            .qbvh
            .refit(0.0, &mut self.workspace, |index| self.aabbs[*index]);
        self.qbvh.rebalance(0.0, &mut self.workspace);
        index
    }

    fn remove_aabb(&mut self, aabb_index: usize) {
        let _aabb = self.aabbs.remove(aabb_index);
        let _removed = self.qbvh.remove(aabb_index);
        let _count = self
            .qbvh
            .refit(0.0, &mut self.workspace, |index| self.aabbs[*index]);
        self.qbvh.rebalance(0.0, &mut self.workspace);
    }

    fn check_topology(&self) {
        self.qbvh.check_topology(true, |index| {
            *self
                .aabbs
                .get(*index)
                .expect(&format!("invalid index {}", index))
        });
    }

    fn print_tree(&self) {
        let mut string = vec![];
        let _ = ptree::write_tree(&QbvhTreeIterator::new(&self.qbvh), &mut string);
        println!("{}", String::from_utf8(string).unwrap());
    }
}

fn generate_random_aabb(rng: &mut StdRng) -> Aabb {
    let min_x = rng.gen_range(-100..100);
    let min_y = rng.gen_range(-100..100);
    #[cfg(feature = "dim3")]
    let min_z = rng.gen_range(-100..100);

    let max_x = rng.gen_range(min_x..(min_x + 50));
    let max_y = rng.gen_range(min_y..(min_y + 50));
    #[cfg(feature = "dim3")]
    let max_z = rng.gen_range(min_z..(min_z + 50));

    #[cfg(feature = "dim3")]
    {
        let mins = Point::new(min_x as Real, min_y as Real, min_z as Real);
        let maxs = Point::new(max_x as Real, max_y as Real, max_z as Real);
        Aabb::new(mins, maxs)
    }
    #[cfg(feature = "dim2")]
    {
        let mins = Point::new(min_x as Real, min_y as Real);
        let maxs = Point::new(max_x as Real, max_y as Real);
        Aabb::new(mins, maxs)
    }
}

impl<'q> ptree::TreeItem for QbvhTreeIterator<'q> {
    type Child = Self;

    fn write_self<W: std::io::Write>(
        &self,
        f: &mut W,
        _style: &ptree::Style,
    ) -> std::io::Result<()> {
        #[cfg(feature = "dim3")]
        return writeln!(
            f,
            "{}{} [{:>4.0}, {:>4.0}, {:>4.0}] -> [{:>4.0}, {:>4.0}, {:>4.0}]",
            if self.child_of_leaf { "*" } else { " " },
            self.node,
            self.aabb.mins.x,
            self.aabb.mins.y,
            self.aabb.mins.z,
            self.aabb.maxs.x,
            self.aabb.maxs.y,
            self.aabb.maxs.z,
        );
        #[cfg(feature = "dim2")]
        return writeln!(
            f,
            "{}{} [{:>4.0}, {:>4.0}] -> [{:>4.0}, {:>4.0}]",
            if self.child_of_leaf { "*" } else { " " },
            self.node,
            self.aabb.mins.x,
            self.aabb.mins.y,
            self.aabb.maxs.x,
            self.aabb.maxs.y,
        );
    }

    fn children(&self) -> Cow<[Self::Child]> {
        self.get_children_node()
    }
}

#[derive(Clone)]
struct QbvhTreeIterator<'q> {
    qbvh: &'q GenericQbvh<usize, DefaultStorage>,
    node: u32,
    aabb: Aabb,
    child_of_leaf: bool,
}

impl<'q> QbvhTreeIterator<'q> {
    fn new(qbvh: &'q GenericQbvh<usize, DefaultStorage>) -> Self {
        let aabb = *qbvh.root_aabb();
        let child_of_leaf = qbvh.raw_nodes().is_empty();

        Self {
            qbvh,
            node: 0,
            aabb,
            child_of_leaf,
        }
    }

    fn get_children_node(&self) -> Cow<[Self]> {
        if self.child_of_leaf {
            return Cow::Borrowed(&[]);
        }

        if let Some(node) = self.qbvh.raw_nodes().get(self.node as usize) {
            let self_leaf = node.is_leaf();

            node.children
                .iter()
                .enumerate()
                .filter_map(|(child_index, &child)| {
                    if child != u32::MAX {
                        let simd_aabb = self.qbvh.raw_nodes()[self.node as usize].simd_aabb;
                        let aabb = simd_aabb.extract(child_index);
                        Some(Self {
                            qbvh: self.qbvh,
                            node: child,
                            aabb,
                            child_of_leaf: self_leaf,
                        })
                    } else {
                        None
                    }
                })
                .collect()
        } else {
            // ?
            panic!("Unknown node index")
        }
    }
}
