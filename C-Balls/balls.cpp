#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>

using namespace sf;
using namespace std;

class Particle {
private:
    float radius;
    int* color;
    Vector2f pos;
    Vector2f last_pos;
    Vector2f accel;

    CircleShape shape;

public:
    Particle() {
        this->radius = 10;
        int temp[3] = { 255, 255, 255 };
        this->color = temp;
        this->pos = Vector2f(10.f, 10.f);
        this->last_pos = Vector2f(10.f, 10.f);
        this->accel = Vector2f(0.f, 0.f);

        this->shape.setFillColor(Color((uint32_t)color[0], (uint32_t)color[1], (uint32_t)color[2]));
        this->shape.setPosition(pos);
        this->shape.setRadius(radius);
        this->shape.setOrigin(Vector2f(radius, radius));
    }

    Particle(float radius, int* color, Vector2f pos, Vector2f accel, Vector2f velocity) {
        this->radius = radius;
        this->color = color;
        this->pos = pos;
        this->last_pos = Vector2f(pos.x - velocity.x, pos.y - velocity.y);
        this->accel = accel;

        this->shape.setFillColor(Color(color[0], color[1], color[2]));
        this->shape.setPosition(pos);
        this->shape.setRadius(radius);
		this->shape.setOrigin(Vector2f(radius, radius));
    }

	void update_padding(float padding) {
		//cout << "Padding: " << padding << endl;
		this->shape.setRadius(this->radius + padding);
		/*this->shape.setOrigin(Vector2f(this->radius + padding, this->radius + padding));*/
	}

    CircleShape render() {
        this->shape.setPosition(this->pos);
        //this->speed_to_color();

        return this->shape;
    }

    void update(float dt) {
        Vector2f displacement = Vector2f(this->pos.x - this->last_pos.x, this->pos.y - this->last_pos.y);

        this->last_pos = Vector2f(this->pos.x, this->pos.y);
        float result = (dt * dt);

        this->pos += Vector2f(displacement.x + (this->accel.x * result), displacement.y + (this->accel.y * result));
        this->accel = Vector2f(0.0f, 0.0f);
        capVelocity();
    }

    void speed_to_color() {
        int intensity = min(255, int((abs(this->getVelocity().lengthSquared()) / 4) * 255));
        float speed = this->getVelocity().lengthSquared();
        float normalized = min(speed / 4.0f, 1.f);
        // 36
        // 167
        int targetR = 255;
		int targetG = 36;
		int targetB = 167;

		int* color = new int[3];
		color[0] = max(int(targetR * normalized), 1);
        color[1] = max(int(targetG * normalized), 1);
        color[2] = max(int(targetB * normalized), 1);
		this->shape.setFillColor(Color(color[0], color[1], color[2]));
    }

    void pressure_to_color(float pressure) {
        int *color = new int[3];
		color[0] = int(255 * pressure);
        color[1] = 0;
		color[2] = int(255 * (1.f - pressure));
		this->shape.setFillColor(Color(color[0], color[1], color[2]));
    }
    

    void accelerate(Vector2f a) {
        this->accel += a;
    }

    void addVelocity(Vector2f v, float dt) {
        this->last_pos -= v * dt;
    }

    void setVelocity(Vector2f v, float dt) {
        this->last_pos = this->pos - (v * dt);
    }

	void setPosition(Vector2f pos) {
		this->pos = pos;
	}

    void capVelocity() {
		Vector2f velocity = this->getVelocity();
		float max_velocity = 5.f; // Set a maximum velocity
		if (velocity.lengthSquared() > max_velocity * max_velocity) {
			velocity = velocity.normalized() * max_velocity;
			this->setVelocity(velocity, 1.f);
		}
    }

    Vector2f getVelocity() const {
        return this->pos - this->last_pos;
    }

    Vector2f getPosition() const {
        return this->pos;
    }

	float getRadius() const {
		return this->radius;
	}
};

class Solver {
private:

	vector<Particle> particles;
	Vector2f gravity = Vector2f(0.f, 1000.f);
    Vector2i window_size;

    int subsets = 8;
    float sub_step_dt = 1.f / (60.f * subsets);
	const float GRID_SIZE = 6.f;

    float distance = 0.f;
    float radius_sum = 0.f;
    Vector2f velo;

public:
	// Constructor
	Solver(vector<Particle> particles, Vector2i window_size) {
		this->particles = particles;
		this->window_size = window_size;
	}

	void update_padding(float padding) {
		for (int i = 0; i < particles.size(); i++) {
			this->particles[i].update_padding(padding);
		}
	}

	void render(RenderWindow &window) {
		for (int i = 0; i < particles.size(); i++) {
			window.draw(this->particles[i].render());
		}
	}

    void updateParticles() {
		for (int i = 0; i < particles.size(); i++) {
			this->particles[i].update(sub_step_dt);
		}
    }

    void applyGravity() {
		for (int i = 0; i < particles.size(); i++) {
			this->particles[i].accelerate(gravity);
		}
    }

	void setParticles(vector<Particle> particles) {
		this->particles = particles;
	}

	void appendParticle(Particle p) {
		this->particles.push_back(p);
	}

    float normalize(float value, float min, float max) {
        return (max - min) > 0 ? (value - min) / (max - min) : 0.1f;
    }

    void applyCircleBoundary() {
    

        for (int i = 0; i < particles.size(); i++) {
            Vector2f pos = this->particles[i].getPosition();
            Vector2f velocity = this->particles[i].getVelocity();
            float radius = this->particles[i].getRadius();
            float circle_radius = 500.f;

            // Center of the window
            Vector2f center(window_size.x / 2.f, window_size.y / 2.f);

            // Vector from center to particle
            Vector2f to_particle = pos - center;
            float dist_sq = to_particle.x * to_particle.x + to_particle.y * to_particle.y;

            if (dist_sq > (circle_radius - radius) * (circle_radius - radius)) {
                float angle = atan2(to_particle.y, to_particle.x);
                Vector2f new_pos = center + Vector2f(
                    cos(angle) * (circle_radius - radius),
                    sin(angle) * (circle_radius - radius)
                );
                this->particles[i].setPosition(new_pos);
                this->particles[i].setVelocity(Vector2f(-velocity.x, -velocity.y), 1.0f);
            }
        }
    }

    void applyBorder() {
        for (int i = 0; i < particles.size(); i++) {
            Vector2f pos = this->particles[i].getPosition();
            Vector2f new_pos = Vector2f(pos.x, pos.y);
            Vector2f velocity = this->particles[i].getVelocity();

            float radius = this->particles[i].getRadius();

            float dampening = 0.75f;

            Vector2f dy = Vector2f(velocity.x * dampening, -velocity.y);
            Vector2f dx = Vector2f(-velocity.x, velocity.y * dampening);

            if (pos.x - radius < 0 || pos.x + radius >= this->window_size.x) {
                if (pos.x - radius < 0) {
                    new_pos.x = radius;
                }
                else {
                    new_pos.x = this->window_size.x - radius;
                }
                this->particles[i].setPosition(new_pos);
                this->particles[i].setVelocity(dx, 1.0f);
            }

            if (pos.y - radius < 0 || pos.y + radius >= this->window_size.y) {
                if (pos.y - radius < 0) {
                    new_pos.y = radius;
                }
                else {
                    new_pos.y = this->window_size.y - radius;
                }
                this->particles[i].setPosition(new_pos);
                this->particles[i].setVelocity(dy, 1.0f);
            }
        }
    }

    void handleCollisions(auto& obj1, auto& obj2) {
        Vector2f normal = this->velo / distance;
        float delta = 0.5 * (this->radius_sum - this->distance);

        obj1.setPosition(obj1.getPosition() + normal * delta);
        obj2.setPosition(obj2.getPosition() - normal * delta);
    }

    bool collide(auto& obj1, auto& obj2) {
        Vector2f pos1 = obj1.getPosition();
        Vector2f pos2 = obj2.getPosition();
        this->velo = pos1 - pos2;

		float dist_squared = velo.x * velo.x + velo.y * velo.y;
        this->radius_sum = obj1.getRadius() + obj2.getRadius();
        float min_dist_squarted = pow(radius_sum, 2);
        
        
        if (dist_squared < min_dist_squarted) {
            this->distance = sqrt(dist_squared);  
            return true;
        }
        return false;

    }

    void find_collisions() {
        unordered_map<int, vector<Particle*>> grid;
		//vector<float> pressures(particles.size(), 0.f);

		// Populate the grid with particles
        for (auto& particle : particles) {
            int cellX = static_cast<int>(std::floor(particle.getPosition().x / GRID_SIZE));
            int cellY = static_cast<int>(std::floor(particle.getPosition().y / GRID_SIZE));
            int cellHash = cellX * 73856093 ^ cellY * 19349663; // Hash function for 2D grid
            grid[cellHash].push_back(&particle);
        }

		// Check for collisions for neighboring cells
        for (auto& particle : particles) {
			Vector2f p_pos = particle.getPosition();
			Vector2f p_grid_pos = Vector2f(p_pos.x / GRID_SIZE, p_pos.y / GRID_SIZE);
			for (int i = -1; i < 2; i++) {
				for (int j = -1; j < 2; j++) {
					int cellX = static_cast<int>(p_grid_pos.x) + j;
					int cellY = static_cast<int>(p_grid_pos.y) + i;
					int cellHash = cellX * 73856093 ^ cellY * 19349663; // Hash function for 2D grid
					if (grid.find(cellHash) != grid.end()) {
						std::vector<Particle*>& cellParticles = grid[cellHash];
						// Check for collisions only within the specific neighbor
						for (size_t k = 0; k < cellParticles.size(); ++k) {
							if (cellParticles[k] == &particle) continue; // Skip self-collision
							if (collide(particle, *cellParticles[k])) {
								handleCollisions(particle, *cellParticles[k]);
							}                       
						}
					}
				}
			}

        }
   
    }  

	void mouse_pull(Vector2f pos, int size) {
		for (int i = 0; i < particles.size(); i++) {
            Vector2f p_pos = this->particles[i].getPosition();
			Vector2f dir = pos - p_pos;

			float distance = sqrt(pow(dir.x, 2) + pow(dir.y, 2));
            this->particles[i].accelerate(dir * max(0.f, 10 * (size - distance)));
		}
	}

	void mouse_push(Vector2f pos, int size) {
        for (int i = 0; i < particles.size(); i++) {
            Vector2f p_pos = this->particles[i].getPosition();
            Vector2f dir = pos - p_pos;

            float distance = sqrt(pow(dir.x, 2) + pow(dir.y, 2));
            this->particles[i].accelerate(dir * min(0.f, -10 * (size - distance)));
        }
	}

	void update(float dt, RenderWindow &window) {
        for (int j = 0; j < subsets; j++) {
            applyGravity();
            updateParticles();

            find_collisions();

            applyBorder();
            //applyCircleBoundary();
            //applyPressure();
        }
        render(window);
	}
 };


 class Solver2 {
 private:
     Vector2f* pos;
     Vector2f* last_pos;
     Vector2f* acceleration;
     float* radius;
     Vector2f gravity = Vector2f(0.f, 1000.f);
     Vector2i window_size;
     int ball_count;
     CircleShape shape;

	 float GRID_SIZE = 6.f;

     float distance = 0.f;
     float radius_sum = 0.f;
     Vector2f velo;

     int subsets = 8;
     float sub_step_dt = 1.f / (60.f * subsets);
 public:
     Solver2(Vector2f* pos, Vector2f* last_pos, Vector2f* acceleration, float* radius, int ball_count, Vector2i window_size) {
         this->pos = pos;
         this->last_pos = last_pos;
         this->acceleration = acceleration;
         this->radius = radius;
         this->ball_count = ball_count;
         this->window_size = window_size;

		 this->shape.setRadius(3.f);
         this->shape.setFillColor(Color(0, 0, 255));
         shape.setOrigin(Vector2f(3.f, 3.f));
     }
	 void render(RenderWindow& window) {
		 for (int i = 0; i < ball_count; i++) {
			 shape.setPosition(pos[i]);
			 window.draw(shape);
		 }
	 }
     void applyGravity() {
		 for (int i = 0; i < ball_count; i++) {
			 this->acceleration[i] += gravity;
		 }
		 //cout << "G" << this->acceleration[0].x << ", " << this->acceleration[0].y << endl;
     }
     void updateParticles() {
         for (int i = 0; i < ball_count; i++) {
             Vector2f displacement = Vector2f(pos[i].x - last_pos[i].x, pos[i].y - last_pos[i].y);
             this->last_pos[i] = Vector2f(pos[i].x, pos[i].y);
             float result = (sub_step_dt * sub_step_dt);
             this->pos[i] += Vector2f(displacement.x + (acceleration[i].x * result), displacement.y + (acceleration[i].y * result));
             this->acceleration[i] = Vector2f(0.f, 0.f);
			 capVelocity(i);
         }
		 /*cout << "Pos: " << pos[0].x << ", " << pos[0].y << endl;
		 cout << "Last Pos: " << last_pos[0].x << ", " << last_pos[0].y << endl;*/
     }
     void setVelocity(Vector2f v, float dt, int idx) {
        this->last_pos[idx] = this->pos[idx] - (v * dt);
     }
     void capVelocity(int i) {
         Vector2f velocity = this->pos[i] - this->last_pos[i];
         float max_velocity = 5.f; // Set a maximum velocity
         if (velocity.lengthSquared() > max_velocity * max_velocity) {
             velocity = velocity.normalized() * max_velocity;
             setVelocity(velocity, 1.f, i);
         }
     }

     void applyBorder() {
         for (int i = 0; i < ball_count; i++) {
             Vector2f pos = Vector2f(this->pos[i].x , this->pos[i].y);
             Vector2f new_pos = Vector2f(pos.x, pos.y);
             Vector2f velocity = Vector2f(this->pos[i].x - this->last_pos[i].x, this->pos[i].y - this->last_pos[i].y);

             float radius = this->radius[i];

             float dampening = 0.75f;

             Vector2f dy = Vector2f(velocity.x * dampening, -velocity.y);
             Vector2f dx = Vector2f(-velocity.x, velocity.y * dampening);

             if (pos.x - radius < 0 || pos.x + radius >= this->window_size.x) {
                 if (pos.x - radius < 0) {
                     new_pos.x = radius;
                 }
                 else {
                     new_pos.x = this->window_size.x - radius;
                 }
                 this->pos[i] = new_pos;
                 //cout << "Last Pos Before (SetV)" << last_pos[i].x << ", " << last_pos[i].y << endl;
                 setVelocity(dx, 1.0f, i);
				 //cout << "Last Pos After (SetV)" << last_pos[i].x << ", " << last_pos[i].y << endl;
             }

             if (pos.y - radius < 0 || pos.y + radius >= this->window_size.y) {
                 if (pos.y - radius < 0) {
                     new_pos.y = radius;
                 }
                 else {
                     new_pos.y = this->window_size.y - radius;
                 }
                 this->pos[i] = new_pos;
                 //cout << "Last Pos Before (SetV)" << last_pos[i].x << ", " << last_pos[i].y << endl;
                 setVelocity(dy, 1.0f, i);
                 //cout << "Last Pos After (SetV)" << last_pos[i].x << ", " << last_pos[i].y << endl;
             }
         }
     }
     void handleCollisions(int idx_1, int idx_2) {
         Vector2f normal = this->velo / distance;
         float delta = 0.5 * (this->radius_sum - this->distance);
		 //cout << "Before" << endl;
		 //cout << "Pos 1: " << this->pos[idx_1].x << ", " << this->pos[idx_1].y << endl;
		 //cout << "Pos 2: " << this->pos[idx_2].x << ", " << this->pos[idx_2].y << endl;
         this->pos[idx_1] += (normal * delta);
         this->pos[idx_2] -= (normal * delta);
		 //cout << "After" << endl;
		 //cout << "Pos 1: " << this->pos[idx_1].x << ", " << this->pos[idx_1].y << endl;
		 //cout << "Pos 2: " << this->pos[idx_2].x << ", " << this->pos[idx_2].y << endl;
     }

     bool collide(int idx_1, int idx_2) {
         Vector2f pos1 = pos[idx_1];
         Vector2f pos2 = pos[idx_2];
         this->velo = pos[idx_1] - pos[idx_2];
		 //cout << "Velo: " << velo.x << ", " << velo.y << endl;
         float dist_squared = (velo.x * velo.x) + (velo.y * velo.y);
		 //cout << "Dist Squared: " << dist_squared << endl;
         this->radius_sum = radius[idx_1] + radius[idx_2];
         float min_dist_squarted = pow(radius_sum, 2);


         if (dist_squared < min_dist_squarted) {
             this->distance = sqrt(dist_squared);
             return true;
         }
         return false;
     }

     void find_collisions() {
         unordered_map<int, vector<int>> grid;

         // Populate the grid with particles
         for (int i = 0; i < ball_count; i++) {
             int cellX = static_cast<int>(std::floor(pos[i].x / GRID_SIZE));
             int cellY = static_cast<int>(std::floor(pos[i].y / GRID_SIZE));
             int cellHash = cellX * 73856093 ^ cellY * 19349663; // Hash function for 2D grid
             //cout << "CellHash: " << cellHash << " I: " << i << endl;
             grid[cellHash].push_back(i);
         }

         // Check for collisions for neighboring cells
         for (int o = 0; o < ball_count; o++) {
             Vector2f p_pos = Vector2f(pos[o].x, pos[o].y);
             Vector2f p_grid_pos = Vector2f(p_pos.x / GRID_SIZE, p_pos.y / GRID_SIZE);
             for (int i = -1; i < 2; i++) {
                 for (int j = -1; j < 2; j++) {
                     int cellX = static_cast<int>(p_grid_pos.x) + j;
                     int cellY = static_cast<int>(p_grid_pos.y) + i;
                     int cellHash = cellX * 73856093 ^ cellY * 19349663; // Hash function for 2D grid
                     if (grid.find(cellHash) != grid.end()) {
                         std::vector<int> cellParticles = grid[cellHash];
                         // Check for collisions only within the specific neighbor
                         for (size_t k = 0; k < cellParticles.size(); k++) {
                             if (cellParticles[k] == o) {
								 /*cout << "Skipping Self Collision" << o << " " << cellParticles[k] << endl;*/
                                 continue; // Skip self-collision
                             }
                             if (collide(o, cellParticles[k])) {
                                 //cout << "Hit" << " o: " << o << " cellK: " << cellParticles[k] << endl;
                                 handleCollisions(o, cellParticles[k]);
                             }
                         }
                     }
                 }
             }

         }

     }


	 void update(float dt, RenderWindow& window) {
         for (int i = 0; i < subsets; i++) {
			 applyGravity();
             updateParticles();

             find_collisions();

			 applyBorder();

         }

		 
		 render(window);
	 }
 };



int main()
{
    sf::ContextSettings settings;
    settings.antiAliasingLevel = 8;
    
    sf::RenderWindow window(sf::VideoMode({ 1000, 1000 }), "Ball Sim", sf::Style::Default, sf::State::Windowed, settings);
    window.setFramerateLimit(60);

    Clock clock;

    vector<Particle> balls;
    int ball_count = 0;

	int* color = new int[3];
	color[0] = 0;
    color[1] = 0;
    color[2] = 255;
    for (int i = 0; i < 10000; i++) {
        //balls.push_back(Particle(3.f, color, Vector2f(10, 10 + 20*i), Vector2f()));
        balls.push_back(Particle(3.f, color, Vector2f(rand() % window.getSize().x, rand() % window.getSize().y), Vector2f(), Vector2f()));
		//balls[i].setVelocity(Vector2f(rand() % 2, rand() % 2), 1.f);
		ball_count++;
    }

	Vector2f* pos = new Vector2f[ball_count];  
	Vector2f* last_pos = new Vector2f[ball_count];
	Vector2f* acceleration = new Vector2f[ball_count];
	float* radius = new float[ball_count];

	for (int i = 0; i < ball_count; i++) {
		pos[i] = balls[i].getPosition();
		last_pos[i] = Vector2f(pos[i].x - balls[i].getVelocity().x, pos[i].y - balls[i].getVelocity().y);
		acceleration[i] = Vector2f(0.f, 0.f);
		radius[i] = balls[i].getRadius();
	}

	Solver2 solver2(pos, last_pos, acceleration, radius, ball_count, (Vector2i)window.getSize());

    Solver solver(balls, (Vector2i)window.getSize());
	
    int frameCount = 0;
	bool stop_spawning = true;

    int size = 120;
    float padding = 0.f;
    bool ctrlHeld = false;
    // run the program as long as the window is open
    while (window.isOpen())
    {

        // check all the window's events that were triggered since the last iteration of the loop
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>())
                window.close();

            if (event->is<sf::Event::MouseWheelScrolled>()) {
                const auto& scroll = event->getIf<Event::MouseWheelScrolled>();
                //Event::KeyPressed key;
                if (ctrlHeld) {
                    if (padding > 0 && scroll->delta < 0) {
                        padding -= 0.5f;
                    }
                    else if (padding < 4.f && scroll->delta > 0) {
                        padding += 0.5f;
                    }
                    //cout << "Padding: " << padding << endl;
                    solver.update_padding(padding);

                }
                else {
                    if (scroll->delta < 0) {
                        if (size > 80) {
                            size -= 5;
                        }
                    }
                    else {
                        if (size < 180) {
                            size += 5;
                        }
                    }
                }
            }
            //cout << "Size: " << size << endl;


            if (event->is<sf::Event::KeyPressed>()) {
                const auto& key = event->getIf<Event::KeyPressed>();
                if (key->code == sf::Keyboard::Key::R) {
                    balls.clear();
                    ball_count = 0;
                    for (int i = 0; i < 8000; i++) {
                        balls.push_back(Particle(3.f, color, Vector2f(rand() % window.getSize().x, rand() % window.getSize().y), Vector2f(), Vector2f()));
                        //balls[i].setVelocity(Vector2f(rand() % 2, rand() % 2), 1.f);
                        ball_count++;
                    }
                    solver.setParticles(balls);
                    solver.update_padding(padding);
                    //stop_spawning = false;
                }
                if (key->code == sf::Keyboard::Key::Escape) {
                    window.close();
                }

				if (key->code == sf::Keyboard::Key::LControl) {
					ctrlHeld = true;
				}
            }

            if (event->is<Event::KeyReleased>()) {
				const auto& key = event->getIf<Event::KeyReleased>();
				if (key->code == sf::Keyboard::Key::LControl) {
                    ctrlHeld = false;
				}
            }
        
        }

        float currentTime = clock.restart().asSeconds();
        float fps = 1.0f / currentTime;
        window.setTitle("Ball Sim - P: " + to_string(ball_count) + " - FPS: " + to_string(fps));
        
        if (fps < 30) {
            stop_spawning = true;
        }
        if (!stop_spawning) {
            solver.appendParticle(Particle(3.f, color, Vector2f(10, 10), Vector2f(), Vector2f(1, 0)));
            ball_count++;
        }
		
		if (Mouse::isButtonPressed(Mouse::Button::Left)) {
			Vector2f m_pos = (Vector2f) Mouse::getPosition(window);
            solver.mouse_pull(m_pos, size);
		}

        if (Mouse::isButtonPressed(Mouse::Button::Right)) {
            Vector2f m_pos = (Vector2f)Mouse::getPosition(window);
            solver.mouse_push(m_pos, size);
        }        
		//cout << "FPS: " << fps << endl;

		// clear the window with black color    
        window.clear(sf::Color::Black);

		//solver2.update(1.f / 60.f, window);
		solver.update(1.f / 60.f, window);

        window.display();
    }
    return 0;
}