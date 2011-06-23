/* Chapel program solving the N-Body problem. */

config const G = 6.6748e-11;
config const numBodies = 2;
config const numDimensions = 2;
config const inputFileName = "default.in";
config const timeStep = 20.0;
config const timeFinal = 2e4;


record OneBody {
    var mass: real;
    var position: [1..numDimensions] real;
    var velocity: [1..numDimensions] real;
}

proc distance(x: [],y: []) {
    return sqrt(+ reduce (x - y)**2);
}

class PhysicalProblem {
    // abstract class
    var time: real;
    proc get_y() {};
    proc update(y:real) {};
}

class NBodyProblem : PhysicalProblem {
    var points: [1..numBodies] OneBody;
    const stride = numBodies * numDimensions;

    proc readFromFile(filename) {
        var infile = new file(filename, FileAccessMode.read);
        infile.open();

        for p in points {
            infile.read(p.mass);
            for j in 1..numDimensions do 
                infile.read(p.position[j]);
            for j in 1..numDimensions do
                infile.read(p.velocity[j]);
        }
        infile.close(); 
    }

    proc printParameters(){
         for p in points {
            writeln("mass: ", p.mass);
            writeln(" position ", p.position);
            writeln(" velocity ", p.velocity);
        }
    }
    
    proc printActualStatus() {
        writeln(time, " ", get_y());
    }

    proc acceleration() {
        var a: [1..stride] real;
        var allBodies = 1..numBodies;

        forall (i,j) in [allBodies, allBodies] { // tensor iteration
            if i != j {
                var distance_vec = points[i].position - points[j].position;
                var r = distance(points[i].position, points[j].position);
                var ind_i = (1..numDimensions) + (i-1)*numDimensions;

                a[ind_i] -= G * points[j].mass * distance_vec / r**3;
            }
        }

        return a;
    }
    
    proc rhs(y, t: real) {
        var rhs: [1..2*stride] real;
        var tmp: [1..2*stride] real;

        var lower = [1..stride];
        var upper = [stride+1..2*stride];

        tmp = y;
        update(tmp);

        rhs[lower] = tmp[upper];
        rhs[upper] = acceleration();

        return rhs;
    }
    
    proc get_y() {
        var y: [1 .. 2 * stride] real;
        var k = 1;
        for p in points do {
            for d in 1..numDimensions do {
                y[k] = p.position[d];
                y[k + stride] = p.velocity[d];
                k += 1;
            }
        }
        return y;
    }
    
    proc update(y: [] real) {
        var k=1;
        for p in points do {
            for d in 1..numDimensions do {
                p.position[d] = y[k];
                p.velocity[d] = y[k + stride];
                k += 1;
            }
        }
    }
}


class RungeKutta4 {
    var problem : NBodyProblem;
    proc advance(dt: real) {
        var y = problem.get_y();
        var t = problem.time;

        var k1 = problem.rhs(y, t);
        var k2 = problem.rhs(y + 0.5*dt*k1, t+0.5*dt);
        var k3 = problem.rhs(y + 0.5*dt*k2, t+0.5*dt);
        var k4 = problem.rhs(y + dt*k3, t+dt);

        var y_new = y + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);

        problem.update(y_new); 
        problem.time += dt;
    }
}


proc main {
    var problem = new NBodyProblem(); 

    problem.readFromFile(inputFileName);

    var step_engine = new RungeKutta4(problem);
    
    while problem.time <= timeFinal do {
        problem.printActualStatus();
        step_engine.advance(timeStep);
    }
}

