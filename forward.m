function results = forward(N, num_particles, init_fn, choice_fn, update_fn)

    choices = [];
    for i=1:num_particles
        particles(i) = init_fn();
        [w(i) c(i)] = choice_fn(1, particles(i));
    end
    choices = [choices; mean(c)];

    particles = resample_particles(particles, w);

    for n=2:N
        for i=1:num_particles
            particles(i) = update_fn(n-1, particles(i));
            [w(i) c(i)] = choice_fn(n, particles(i));
        end
        choices = [choices; mean(c)];
    
        particles = resample_particles(particles, w);
    end
    
    for i=1:num_particles
        particles(i) = update_fn(n, particles(i)); % last update of posterior
    end

    results.ww_n = particles(1).w;
    results.choices = choices;
    results.P_n = mean(cat(1,particles.sample), 1); % for model_test.m
    results.sample = mean(cat(1,particles.sample), 1); 
    results.samples = mean(cat(3, particles.samples), 3);
    results.Posterior = mean(cat(3, particles.Posterior), 3);


