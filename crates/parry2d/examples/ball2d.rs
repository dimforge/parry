use parry2d::shape::Ball;

fn main() {
    let ball = Ball::new(1.0);
    assert!(ball.radius == 1.0);
}
