using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class particle
{
	public Transform trans;
	public float mass;
	Vector3 acceleration;
	public Vector3 velocity;
	public float d0;

	float viscosity = 2.5f;

	public Vector3 position
	{
		get
		{
			return trans.position;
		}
		set
		{
			trans.position = value;
		}
	}

	public particle(Transform _trans, float _mass, float _d0)
	{
		trans = _trans;
		mass = _mass;
		d0 = _d0;
	}

	public void ApplyForce(Vector3 force)
	{
		acceleration = force / mass;
		velocity = velocity + Time.deltaTime * acceleration;
		position = position + Time.deltaTime * velocity;
	}

	public void PositiveVelocityY()
	{
		velocity = new Vector3(velocity.x / 1.1f, Mathf.Abs(velocity.y) / (3 * viscosity), velocity.z);
	}

	public void PositiveVelocityX()
	{
		velocity = new Vector3(Mathf.Abs(velocity.x) / viscosity, velocity.y / (0.9f * viscosity), velocity.z);
	}

	public void NegativeVelocityX()
	{
		velocity = new Vector3(-Mathf.Abs(velocity.x) / viscosity, velocity.y / (0.9f * viscosity), velocity.z);
	}
}

public class CreateParticles : MonoBehaviour
{
	[SerializeField] float gravity = 2f;
	[SerializeField] float rayon_particule = 0.1f;
	[SerializeField] int particle_number = 100;
	[SerializeField] GameObject particle_demo;
	[SerializeField] GameObject particle_2_demo;

	float width_create;
	float height_create;

	float boundX_screen;
	float boundY_screen;
	float rayon_particule_real = 0.1f;

	[SerializeField] float k0 = 2.0f;
	[SerializeField] float k1 = 0.5f;
	[SerializeField] float d0_1 = 3f;
	[SerializeField] float d0_2 = 0.8f;

	float step;

	List<particle> particles = new List<particle> ();

	Dictionary<int, List<particle>> grid = new Dictionary<int, List<particle>> ();

	// Use this for initialization
	void Start ()
	{
		step = rayon_particule / 1.5f;
		width_create = 1.5f;
		height_create = 1;

		boundX_screen = 2.5f;
		boundY_screen = 5;

		for (int k = 0; k < particle_number; k++)
		{
			SpawnParticle ();
		}
	}

	// Update is called once per frame
	void Update ()
	{
		//SpawnParticle ();

		for (int k = 0; k < particles.Count; k++)
		{
			particle part = particles [k];

			// apply gravity: vi = vi + dt*g
			part.velocity = part.velocity + Vector3.down * Time.deltaTime * gravity;

			// save previous position: xi_prev = xi
			Vector3 pos_old = part.position;
			// advance to predicted position: xi = xi + dt * vi
			part.position = part.position + Time.deltaTime * part.velocity;

			part = DoubleDensityRelaxation(part);

			part = ManageBorders(part);

			// use previous position to compute next velocity: vi = (xi - xi_prev) / dt
			part.velocity = (part.position - pos_old) / Time.deltaTime;
		}
		grid.Clear ();
		for (int k = 0; k < particles.Count; k++)
		{
			UpdateGrid (particles [k]);
		}
	}

	void SpawnParticle()
	{
		if (particles.Count < particle_number)
		{
			GameObject b;
			b = Instantiate (
				particle_demo, 
				new Vector3 (
					transform.position.x + Random.Range (-width_create + 0.2f, width_create - 0.2f), 
					transform.position.y + Random.Range (-height_create, height_create),
					0),
				Quaternion.identity, 
				transform);
			particle part = new particle (b.transform, 0.2f, d0_1);
			part.ApplyForce(new Vector3(Random.Range(-5, 5), 0, 0));
			particles.Add (part);

			UpdateGrid (part);

			b = Instantiate (
				particle_2_demo, 
				new Vector3 (
					transform.position.x + Random.Range (-width_create + 0.2f, width_create - 0.2f), 
					transform.position.y + Random.Range (-height_create, height_create),
					0),
				Quaternion.identity, 
				transform);
			part = new particle (b.transform, 0.2f, d0_2);
			part.ApplyForce(new Vector3(Random.Range(-5, 5), 0, 0));
			particles.Add (part);

			UpdateGrid (part);
		}
	}

	particle ManageBorders(particle part)
	{
		// floor (down)
		if (part.position.y - rayon_particule_real < -boundY_screen)
		{
			float deltaPos = part.position.y - rayon_particule_real - (-boundY_screen);

			part.position = new Vector3 (part.position.x, -boundY_screen, part.position.z);

			part.position += (20f * Vector3.down * deltaPos) * Time.deltaTime * Time.deltaTime;
		}
		// left wall
		if (part.position.x - rayon_particule_real < -boundX_screen)
		{
			float deltaPos = part.position.x - rayon_particule_real + boundX_screen;

			part.position = new Vector3(-boundX_screen * 0.99f, part.position.y, part.position.z);

			part.position += (20f * Vector3.left * deltaPos) * Time.deltaTime * Time.deltaTime;
		}
		// right wall
		else if (part.position.x + rayon_particule_real > boundX_screen)
		{
			float deltaPos = part.position.x + rayon_particule_real - boundX_screen;

			part.position = new Vector3(boundX_screen * 0.99f, part.position.y, part.position.z);

			part.position += (20f * Vector3.left * deltaPos) * Time.deltaTime * Time.deltaTime;
		}

		return part;
	}

	particle DoubleDensityRelaxation(particle part)
	{
		Vector3 moveall = Vector3.zero;

		// get the indices of the virtual box
		int j = Mathf.FloorToInt (part.position.x / step);
		int i = Mathf.FloorToInt (part.position.y / step);

		float density = 0;
		float density2 = 2;
		Vector2 positionStep;

		// compute density and near-density
		for (int m = -1; m <= +1; m++)
		{
			for (int n = -1; n <= +1; n++)
			{
				positionStep = new Vector2 (i + n, j + m);
				// density_i
				density += densityList (positionStep, part, 2);
				// density_i_near
				density2 += densityList (positionStep, part, 3);
			}
		}

		// particle has some neighbouring particles
		if (density > 0)
		{
			// compute pressure and near-pressure
			for (int m = -1; m <= +1; m++)
			{
				for (int n = -1; n <= +1; n++)
				{
					positionStep = new Vector2 (i + n, j + m);
					moveall += moveList (positionStep, part, density, density2);
				}
			}
		}

		// xi = xi + dx
		part.position += moveall;

		return part;
	}

	void UpdateGrid(particle part)
	{
		int j_pos = Mathf.FloorToInt (part.position.x / step);
		int i_pos = Mathf.FloorToInt (part.position.y / step);
		int hash = new Vector2 (i_pos, j_pos).GetHashCode ();

		if (!grid.ContainsKey (hash))
		{
			grid.Add (hash, new List<particle> ());
		}
		grid [hash].Add (part);
	}

	float densityList (Vector2 positionStep, particle part, int pow)
	{
		List<particle> parts;
		int hash = positionStep.GetHashCode ();
		if (grid.ContainsKey (hash))
		{
			parts = grid [hash];
		}
		else
		{
			return 0;
		}
		float density = 0;
		foreach (particle part2 in parts)
		{
			// distance h
			float dist = Vector3.Distance (part2.position, part.position);
			// di = sum{j in Neighbours(i)} ((1 - rij / h) ^ pow)
			if (dist < rayon_particule && dist > 0)
			{
				density += Mathf.Pow((1 - dist / rayon_particule), pow);
			}
		}
		return density;
	}

	Vector3 moveList (Vector2 positionStep, particle part, float dk, float d1k)
	{
		List<particle> parts;
		int hash = positionStep.GetHashCode ();
		if (grid.ContainsKey (hash))
		{
			parts = grid [hash];
		}
		else
		{
			return Vector3.zero;
		}
		Vector3 move = Vector3.zero;
		foreach (particle part2 in parts)
		{
			float dist = Vector3.Distance (part2.position, part.position);
			if (dist < rayon_particule && dist > 0)
			{
				Vector3 rij = (part2.position - part.position).normalized;
				float q = dist / rayon_particule;

				// compute pressure and near-pressure
				float P = k0 * (dk - part.d0);
				float P2 = k0 * (dk - part2.d0);
				float P_near = k1 * (d1k);

				// apply displacements
				//D = dt^2 * (P(1-q)+Pnear(1-q)^2) * rij
				Vector3 D = Time.deltaTime * Time.deltaTime * (P * (1 - q) + P_near * (1 - q) * (1 - q)) * rij;
				Vector3 D2 = Time.deltaTime * Time.deltaTime * (P2 * (1 - q) + P_near * (1 - q) * (1 - q)) * rij;

				//dx =  dx - D/2
				move -= D / 2;
				//xj = xj + D/2
				part2.position += D2 / 2;

				// change velocity of other particle also
				Vector3 V_other = Time.deltaTime * (P2 * (1 - q) + P_near * (1 - q) * (1 - q)) * rij;
				part2.velocity += V_other / 2;
			}
		}
		return move;
	}
}
