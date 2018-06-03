extern crate rand;

#[derive(Clone)]
struct Daisy {
    albedo: f64,
    phenotype_volatility: f64,
}

impl Daisy {

    fn black() -> Daisy {
        Daisy{
            albedo: 0.25,
            phenotype_volatility: 0.05,
        }
    }

    fn white() -> Daisy {
        Daisy{
            albedo: 0.75,
            phenotype_volatility: 0.05,
        }
    }

    fn reproduce_prob(&self, temperature: f64) -> f64 {
        let threshold = 17.5;
        let optimal_temperature = 295.5;
        if (temperature - optimal_temperature).abs() < threshold {
            1.0 - ((temperature - optimal_temperature) / threshold).powi(2)
        } else {
            0.0
        }
    }

    fn offspring(&self) -> Daisy {
        let new_albedo = self.albedo +
            self.phenotype_volatility * (rand::random::<f64>() - 0.5);

        Daisy{
            albedo: if new_albedo < 0.0 { 0.0 }
                    else if new_albedo > 1.0 { 1.0 }
                    else { new_albedo },
            phenotype_volatility: self.phenotype_volatility,
        }
    }

}

#[derive(Clone)]
struct World {
    dim: (usize, usize),
    daisies: Vec<Option<Daisy>>,
    death_rate: f64,
}

impl World {

    fn size(&self) -> usize {
        self.dim.0 * self.dim.1
    }

    fn iter(&self) -> WorldIterator {
        WorldIterator{
            pos: 0,
            world: &self,
        }
    }

    fn at(&self, i: usize) -> Option<&Daisy> {
        if (i >= 0) && (i < self.size()) {
            self.daisies[i].as_ref()
        } else {
            None
        }
    }

    fn above(&self, i: usize) -> Option<&Daisy> {
        if i >= self.dim.1 {
            self.at(i - self.dim.1)
        } else {
            None
        }
    }

    fn below(&self, i: usize) -> Option<&Daisy> {
        if i < self.size() - self.dim.1 {
            self.at(i + self.dim.1)
        } else {
            None
        }
    }

    fn left_of(&self, i: usize) -> Option<&Daisy> {
        if i > 0 {
            self.at(i - 1)
        } else {
            None
        }
    }

    fn right_of(&self, i: usize) -> Option<&Daisy> {
        if i < self.size() - 1 {
            self.at(i + 1)
        } else {
            None
        }
    }

    fn albedo(&self) -> f64 {
        self.daisies.iter()
            .map(|daisy_opt| daisy_opt.as_ref()
                 .map(|d| d.albedo)
                 .unwrap_or(0.5))
            .sum::<f64>() / (self.size() as f64)
    }

    fn temperature_field(&self) -> Vec<f64> {
        let q = 0.125; // neightbour diffusivity
        const S: f64 = 917.0; // solar insolation
        const SB: f64 = 5.67e-8;
        const L: f64 = 1.0; // luminosity

        (0..self.size())
            .map(|i| (1.0 - 4.0*q) * self.at(i).map(|d| d.albedo).unwrap_or(0.5) +
                 q*(self.above(i).map(|d| d.albedo).unwrap_or(0.5) +
                    self.right_of(i).map(|d| d.albedo).unwrap_or(0.5) +
                    self.left_of(i).map(|d| d.albedo).unwrap_or(0.5) +
                    self.below(i).map(|d| d.albedo).unwrap_or(0.5)))
            .map(|albedo| (S * L / SB * (1.0 - albedo)).powf(0.25))
            .collect()
    }

    fn new_randomized(dim: (usize, usize), death_rate: f64) -> World {

        fn new_daisy_opt() -> Option<Daisy> {
            let r= rand::random::<f64>();
            if r < 0.1 {
                Some(Daisy::black())
            } else if r < 0.2 {
                Some(Daisy::white())
            } else {
                None
            }
        }

        World{
            dim: dim.clone(),
            daisies: (0..(dim.0*dim.1)).map(|_| new_daisy_opt()).collect(),
            death_rate: death_rate,
        }
    }

    fn iterate(&self) -> World {

        let mut new_world = self.clone();
        let tfield: Vec<f64> = self.temperature_field();

        // reproduce according to local temperature and species
        for i in 0..self.size() {
            let daisy_opt = self.daisies[i].as_ref();
            if daisy_opt.is_some() {
                continue
            }

            let temp = tfield[i];

            // choose a neighbour at random
            let r: f64 = rand::random();
            let neighbour = if r < 0.25 {
                self.above(i)
            } else if r < 0.5 {
                self.right_of(i)
            } else if r < 0.75 {
                self.below(i)
            } else {
                self.left_of(i)
            };

            let prob = match neighbour {
                Some(_) => neighbour.unwrap().reproduce_prob(temp),
                None => 0.0
            };

            if rand::random::<f64>() < prob {
                new_world.daisies[i] = neighbour.map(|d| d.offspring());
            }
        }

        // die according to global death rate
        for i in 0..self.size() {
            if new_world.daisies[i].is_some() && (rand::random::<f64>() < self.death_rate) {
                new_world.daisies[i] = None;
            }
        }
        new_world
    }

    fn print_stats(&self) -> () {
        let count_empty = self.daisies.iter()
            .filter(|d| d.is_none())
            .map(|_| 1_usize)
            .collect::<Vec<_>>()
            .len();
        let planetary_albedo = self.albedo();

        println!("empty cells: {}", count_empty);
        println!("planetary albedo: {}", planetary_albedo);
    }
}

struct WorldIterator<'a> {
    pos: usize,
    world: &'a World,
}

impl<'a> Iterator for WorldIterator<'a> {
    type Item = &'a Daisy;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.world.at(self.pos);
        self.pos += 1;
        result
    }
}



fn main() {

    println!("it compiled!");
    let mut world = World::new_randomized((20, 20), 0.3);

    for timestep in 0..30 {
        world = world.iterate();
        world.print_stats();
    }

}
