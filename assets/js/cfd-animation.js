// Animazione elegante per il titolo homepage con anime.js
class TitleAnimation {
    constructor() {
        this.initAnimation();
    }
    
    initAnimation() {
        // Attendiamo che anime.js sia caricato
        if (typeof anime === 'undefined') {
            setTimeout(() => this.initAnimation(), 100);
            return;
        }
        
        this.setupTitle();
        this.startAnimation();
    }
    
    setupTitle() {
        const titleElement = document.querySelector('.hero-title');
        if (!titleElement) return;
        
        // Dividi il testo in lettere
        const text = titleElement.textContent;
        titleElement.innerHTML = '';
        
        // Crea span per ogni lettera
        for (let i = 0; i < text.length; i++) {
            const span = document.createElement('span');
            span.textContent = text[i] === ' ' ? '\u00A0' : text[i]; // Non-breaking space
            span.style.display = 'inline-block';
            span.style.opacity = '0';
            span.style.transform = 'translateY(50px) rotateX(-90deg)';
            titleElement.appendChild(span);
        }
    }
    
    startAnimation() {
        const letters = document.querySelectorAll('.hero-title span');
        if (!letters.length) return;
        
        // Timeline principale
        const tl = anime.timeline({
            easing: 'easeOutExpo',
            duration: 1500,
            autoplay: false
        });
        
        // Animazione delle lettere
        tl.add({
            targets: letters,
            opacity: [0, 1],
            translateY: [50, 0],
            rotateX: [-90, 0],
            scale: [0.5, 1],
            delay: anime.stagger(80, {start: 300}),
            duration: 1200,
            easing: 'easeOutElastic(1, .8)'
        });
        
        // Effetto glow
        tl.add({
            targets: '.hero-title',
            textShadow: [
                '0 0 0px rgba(100, 200, 255, 0)',
                '0 0 20px rgba(100, 200, 255, 0.8), 0 0 40px rgba(100, 200, 255, 0.4)'
            ],
            duration: 800,
            easing: 'easeInOutQuad'
        }, '-=600');
        
        // Animazione del sottotitolo
        tl.add({
            targets: '.hero-subtitle',
            opacity: [0, 1],
            translateY: [30, 0],
            duration: 1000,
            easing: 'easeOutQuart'
        }, '-=1000');
        
        // Avvia l'animazione
        tl.play();
        
        // Animazione hover per le lettere
        letters.forEach(letter => {
            letter.addEventListener('mouseenter', () => {
                anime({
                    targets: letter,
                    scale: 1.2,
                    color: '#64c8ff',
                    duration: 300,
                    easing: 'easeOutQuart'
                });
            });
            
            letter.addEventListener('mouseleave', () => {
                anime({
                    targets: letter,
                    scale: 1,
                    color: '#ffffff',
                    duration: 300,
                    easing: 'easeOutQuart'
                });
            });
        });
        
        // Animazione di loop sottile
        this.startLoopAnimation();
    }
    
    startLoopAnimation() {
        const letters = document.querySelectorAll('.hero-title span');
        
        setInterval(() => {
            const randomLetter = letters[Math.floor(Math.random() * letters.length)];
            
            anime({
                targets: randomLetter,
                scale: [1, 1.1, 1],
                rotateZ: [0, 5, 0],
                color: ['#ffffff', '#64c8ff', '#ffffff'],
                duration: 1000,
                easing: 'easeInOutSine'
            });
        }, 3000);
    }
    
    initFields() {
        // Inizializza i campi computazionali
        for (let y = 0; y < this.rows; y++) {
            this.velocityField[y] = [];
            this.pressureField[y] = [];
            this.turbulentKineticEnergy[y] = [];
            this.dissipationRate[y] = [];
            
            for (let x = 0; x < this.cols; x++) {
                // Campo di velocità iniziale con profilo parabolico
                const centerY = this.rows / 2;
                const distFromCenter = Math.abs(y - centerY) / centerY;
                const baseVelocity = 1.0 - distFromCenter * distFromCenter;
                
                this.velocityField[y][x] = {
                    u: baseVelocity * 2.0, // Velocità in x
                    v: 0.0, // Velocità in y
                    magnitude: baseVelocity * 2.0
                };
                
                this.pressureField[y][x] = Math.random() * 0.1;
                this.turbulentKineticEnergy[y][x] = this.turbulentIntensity * baseVelocity * baseVelocity;
                this.dissipationRate[y][x] = Math.pow(this.turbulentKineticEnergy[y][x], 1.5) / this.lengthScale;
            }
        }
    }
    
    initVortices() {
        // Crea vortici coerenti che simulano strutture turbolente grandi
        for (let i = 0; i < this.numVortices; i++) {
            this.vortices.push({
                x: Math.random() * this.width,
                y: Math.random() * this.height,
                strength: (Math.random() - 0.5) * 4.0,
                radius: 30 + Math.random() * 40,
                age: 0,
                maxAge: 500 + Math.random() * 300,
                rotationSpeed: (Math.random() - 0.5) * 0.05,
                vx: (Math.random() - 0.5) * 0.5,
                vy: (Math.random() - 0.5) * 0.3,
                phase: Math.random() * Math.PI * 2
            });
        }
    }
    
    initEddies() {
        // Crea piccoli vortici (eddies) per simulare turbolenza a scala piccola
        for (let i = 0; i < this.numEddies; i++) {
            this.eddies.push({
                x: Math.random() * this.width,
                y: Math.random() * this.height,
                strength: (Math.random() - 0.5) * 2.0,
                radius: 8 + Math.random() * 15,
                age: 0,
                maxAge: 100 + Math.random() * 100,
                rotationSpeed: (Math.random() - 0.5) * 0.1,
                vx: (Math.random() - 0.5) * 1.0,
                vy: (Math.random() - 0.5) * 0.8
            });
        }
    }
    
    
    updateVortexDynamics() {
        // Aggiorna vortici principali (Large Eddy Simulation style)
        this.vortices.forEach(vortex => {
            vortex.x += vortex.vx;
            vortex.y += vortex.vy;
            vortex.age++;
            vortex.phase += vortex.rotationSpeed;
            
            // Oscillazione della forza
            vortex.strength *= (0.999 + 0.001 * Math.sin(vortex.phase));
            
            // Wrapping dei bordi
            if (vortex.x < -vortex.radius) vortex.x = this.width + vortex.radius;
            if (vortex.x > this.width + vortex.radius) vortex.x = -vortex.radius;
            if (vortex.y < -vortex.radius) vortex.y = this.height + vortex.radius;
            if (vortex.y > this.height + vortex.radius) vortex.y = -vortex.radius;
            
            // Rigenerazione del vortice
            if (vortex.age > vortex.maxAge) {
                vortex.x = Math.random() * this.width;
                vortex.y = Math.random() * this.height;
                vortex.strength = (Math.random() - 0.5) * 4.0;
                vortex.age = 0;
                vortex.maxAge = 500 + Math.random() * 300;
                vortex.vx = (Math.random() - 0.5) * 0.5;
                vortex.vy = (Math.random() - 0.5) * 0.3;
            }
        });
        
        // Aggiorna eddies (piccola scala)
        this.eddies.forEach(eddy => {
            eddy.x += eddy.vx;
            eddy.y += eddy.vy;
            eddy.age++;
            
            // Decadimento naturale
            eddy.strength *= 0.998;
            eddy.radius *= 1.001;
            
            // Wrapping
            if (eddy.x < 0) eddy.x = this.width;
            if (eddy.x > this.width) eddy.x = 0;
            if (eddy.y < 0) eddy.y = this.height;
            if (eddy.y > this.height) eddy.y = 0;
            
            // Rigenerazione
            if (eddy.age > eddy.maxAge || eddy.radius > 50) {
                eddy.x = Math.random() * this.width;
                eddy.y = Math.random() * this.height;
                eddy.strength = (Math.random() - 0.5) * 2.0;
                eddy.radius = 8 + Math.random() * 15;
                eddy.age = 0;
                eddy.vx = (Math.random() - 0.5) * 1.0;
                eddy.vy = (Math.random() - 0.5) * 0.8;
            }
        });
    }
    
    updateVelocityField() {
        // Aggiorna il campo di velocità basato su vortici e turbolenza
        for (let y = 0; y < this.rows; y++) {
            for (let x = 0; x < this.cols; x++) {
                const worldX = x * this.cellWidth;
                const worldY = y * this.cellHeight;
                
                let totalU = this.velocityField[y][x].u * 0.95;
                let totalV = this.velocityField[y][x].v * 0.95;
                
                // Influenza dei vortici principali
                this.vortices.forEach(vortex => {
                    const dx = worldX - vortex.x;
                    const dy = worldY - vortex.y;
                    const distance = Math.sqrt(dx * dx + dy * dy);
                    
                    if (distance < vortex.radius * 2) {
                        const influence = Math.exp(-distance / vortex.radius) * vortex.strength;
                        const angle = Math.atan2(dy, dx) + Math.PI / 2;
                        
                        totalU += Math.cos(angle) * influence * 0.02;
                        totalV += Math.sin(angle) * influence * 0.02;
                    }
                });
                
                // Influenza degli eddies
                this.eddies.forEach(eddy => {
                    const dx = worldX - eddy.x;
                    const dy = worldY - eddy.y;
                    const distance = Math.sqrt(dx * dx + dy * dy);
                    
                    if (distance < eddy.radius) {
                        const influence = (1 - distance / eddy.radius) * eddy.strength;
                        const angle = Math.atan2(dy, dx) + Math.PI / 2;
                        
                        totalU += Math.cos(angle) * influence * 0.01;
                        totalV += Math.sin(angle) * influence * 0.01;
                    }
                });
                
                // Aggiunta di rumore turbolento
                const noise = 0.005 * (Math.random() - 0.5);
                totalU += noise;
                totalV += noise;
                
                this.velocityField[y][x].u = totalU;
                this.velocityField[y][x].v = totalV;
                this.velocityField[y][x].magnitude = Math.sqrt(totalU * totalU + totalV * totalV);
                
                // Aggiorna energia cinetica turbolenta
                this.turbulentKineticEnergy[y][x] += 0.01 * (this.velocityField[y][x].magnitude - 1.0);
                this.turbulentKineticEnergy[y][x] = Math.max(0.001, this.turbulentKineticEnergy[y][x] * 0.999);
            }
        }
    }
    
    
    initParticles() {
        for (let i = 0; i < this.numParticles; i++) {
            this.particles.push({
                x: Math.random() * this.width,
                y: Math.random() * this.height,
                vx: 0,
                vy: 0,
                trail: [],
                maxTrailLength: 25 + Math.random() * 15,
                life: Math.random() * 200 + 100,
                maxLife: 300,
                size: 1 + Math.random() * 2,
                color: {
                    r: 100 + Math.random() * 100,
                    g: 150 + Math.random() * 100,
                    b: 200 + Math.random() * 55
                },
                turbulentEnergy: 0
            });
        }
    }
    
    updateParticles() {
        this.particles.forEach(particle => {
            // Interpolazione del campo di velocità
            const col = Math.floor(particle.x / this.cellWidth);
            const row = Math.floor(particle.y / this.cellHeight);
            
            if (col >= 0 && col < this.cols - 1 && row >= 0 && row < this.rows - 1) {
                // Interpolazione bilineare per un movimento più fluido
                const fx = (particle.x / this.cellWidth) - col;
                const fy = (particle.y / this.cellHeight) - row;
                
                const v00 = this.velocityField[row][col];
                const v10 = this.velocityField[row][col + 1];
                const v01 = this.velocityField[row + 1][col];
                const v11 = this.velocityField[row + 1][col + 1];
                
                const u = (1 - fx) * (1 - fy) * v00.u + fx * (1 - fy) * v10.u + 
                         (1 - fx) * fy * v01.u + fx * fy * v11.u;
                const v = (1 - fx) * (1 - fy) * v00.v + fx * (1 - fy) * v10.v + 
                         (1 - fx) * fy * v01.v + fx * fy * v11.v;
                
                // Applica velocità del campo
                particle.vx += u * 0.1;
                particle.vy += v * 0.1;
                
                // Aggiorna energia turbolenta della particella
                const k = (1 - fx) * (1 - fy) * this.turbulentKineticEnergy[row][col] + 
                         fx * (1 - fy) * this.turbulentKineticEnergy[row][col + 1] + 
                         (1 - fx) * fy * this.turbulentKineticEnergy[row + 1][col] + 
                         fx * fy * this.turbulentKineticEnergy[row + 1][col + 1];
                
                particle.turbulentEnergy = k;
            }
            
            // Effetto della turbolenza sulle particelle
            const turbulentForce = Math.sqrt(particle.turbulentEnergy) * 0.05;
            particle.vx += (Math.random() - 0.5) * turbulentForce;
            particle.vy += (Math.random() - 0.5) * turbulentForce;
            
            // Viscosità e dissipazione
            particle.vx *= this.viscosity;
            particle.vy *= this.viscosity;
            
            // Aggiorna posizione
            particle.x += particle.vx;
            particle.y += particle.vy;
            
            // Gestione scia con intensità variabile
            particle.trail.push({ 
                x: particle.x, 
                y: particle.y, 
                intensity: particle.turbulentEnergy 
            });
            if (particle.trail.length > particle.maxTrailLength) {
                particle.trail.shift();
            }
            
            // Wrapping dei bordi
            if (particle.x < 0) particle.x = this.width;
            if (particle.x > this.width) particle.x = 0;
            if (particle.y < 0) particle.y = this.height;
            if (particle.y > this.height) particle.y = 0;
            
            // Rigenerazione delle particelle
            particle.life--;
            if (particle.life <= 0) {
                particle.x = Math.random() * this.width;
                particle.y = Math.random() * this.height;
                particle.vx = 0;
                particle.vy = 0;
                particle.life = particle.maxLife;
                particle.trail = [];
            }
        });
    }
    
    
    drawVortices() {
        // Disegna i vortici principali con effetto glow
        this.vortices.forEach(vortex => {
            const alpha = Math.max(0.1, 1 - vortex.age / vortex.maxAge);
            const radius = vortex.radius * (0.5 + 0.5 * alpha);
            
            // Gradiente radiale per l'effetto vortice
            const gradient = this.ctx.createRadialGradient(
                vortex.x, vortex.y, 0,
                vortex.x, vortex.y, radius
            );
            gradient.addColorStop(0, `rgba(255, 100, 150, ${alpha * 0.8})`);
            gradient.addColorStop(0.7, `rgba(100, 200, 255, ${alpha * 0.4})`);
            gradient.addColorStop(1, 'rgba(100, 200, 255, 0)');
            
            this.ctx.fillStyle = gradient;
            this.ctx.beginPath();
            this.ctx.arc(vortex.x, vortex.y, radius, 0, Math.PI * 2);
            this.ctx.fill();
            
            // Linee di rotazione per evidenziare il vortice
            this.ctx.strokeStyle = `rgba(255, 150, 200, ${alpha * 0.6})`;
            this.ctx.lineWidth = 2;
            for (let i = 0; i < 6; i++) {
                const angle = (i / 6) * Math.PI * 2 + vortex.phase;
                const startRadius = radius * 0.3;
                const endRadius = radius * 0.8;
                
                this.ctx.beginPath();
                this.ctx.arc(vortex.x, vortex.y, startRadius + (endRadius - startRadius) * i / 6, 
                           angle, angle + Math.PI / 4);
                this.ctx.stroke();
            }
        });
    }
    
    drawEddies() {
        // Disegna gli eddies come piccoli vortici turbolenti
        this.eddies.forEach(eddy => {
            const alpha = Math.max(0.05, 1 - eddy.age / eddy.maxAge);
            
            this.ctx.strokeStyle = `rgba(150, 255, 200, ${alpha})`;
            this.ctx.lineWidth = 1;
            this.ctx.beginPath();
            this.ctx.arc(eddy.x, eddy.y, eddy.radius, 0, Math.PI * 2);
            this.ctx.stroke();
            
            // Piccole linee turbolenti
            for (let i = 0; i < 4; i++) {
                const angle = (i / 4) * Math.PI * 2 + eddy.age * 0.1;
                const x = eddy.x + Math.cos(angle) * eddy.radius * 0.7;
                const y = eddy.y + Math.sin(angle) * eddy.radius * 0.7;
                
                this.ctx.beginPath();
                this.ctx.arc(x, y, 2, 0, Math.PI * 2);
                this.ctx.fill();
            }
        });
    }
    
    drawVelocityField() {
        // Visualizza il campo di velocità con colori e intensità
        for (let y = 0; y < this.rows; y += 2) {
            for (let x = 0; x < this.cols; x += 2) {
                const velocity = this.velocityField[y][x];
                const magnitude = Math.min(1, velocity.magnitude / 3);
                
                if (magnitude > 0.1) {
                    const centerX = x * this.cellWidth + this.cellWidth / 2;
                    const centerY = y * this.cellHeight + this.cellHeight / 2;
                    
                    // Colore basato sulla velocità
                    const hue = 240 - magnitude * 120; // Da blu a rosso
                    this.ctx.strokeStyle = `hsla(${hue}, 70%, 60%, ${magnitude * 0.8})`;
                    this.ctx.lineWidth = 1 + magnitude;
                    
                    const angle = Math.atan2(velocity.v, velocity.u);
                    const length = magnitude * 20;
                    const endX = centerX + Math.cos(angle) * length;
                    const endY = centerY + Math.sin(angle) * length;
                    
                    this.ctx.beginPath();
                    this.ctx.moveTo(centerX, centerY);
                    this.ctx.lineTo(endX, endY);
                    this.ctx.stroke();
                    
                    // Freccia
                    const arrowSize = 2 + magnitude * 2;
                    this.ctx.beginPath();
                    this.ctx.moveTo(endX, endY);
                    this.ctx.lineTo(
                        endX - Math.cos(angle - 0.5) * arrowSize,
                        endY - Math.sin(angle - 0.5) * arrowSize
                    );
                    this.ctx.moveTo(endX, endY);
                    this.ctx.lineTo(
                        endX - Math.cos(angle + 0.5) * arrowSize,
                        endY - Math.sin(angle + 0.5) * arrowSize
                    );
                    this.ctx.stroke();
                }
            }
        }
    }
    
    
    drawParticles() {
        this.particles.forEach(particle => {
            const alpha = particle.life / particle.maxLife;
            
            // Disegna la scia con gradiente di intensità
            if (particle.trail.length > 1) {
                this.ctx.lineWidth = 1;
                
                for (let i = 1; i < particle.trail.length; i++) {
                    const trailAlpha = (i / particle.trail.length) * alpha * 0.8;
                    const turbulentBoost = Math.min(1, particle.trail[i].intensity * 5);
                    
                    // Colore basato sull'energia turbolenta
                    const r = particle.color.r + turbulentBoost * 50;
                    const g = particle.color.g + turbulentBoost * 30;
                    const b = particle.color.b;
                    
                    this.ctx.strokeStyle = `rgba(${r}, ${g}, ${b}, ${trailAlpha})`;
                    this.ctx.beginPath();
                    this.ctx.moveTo(particle.trail[i-1].x, particle.trail[i-1].y);
                    this.ctx.lineTo(particle.trail[i].x, particle.trail[i].y);
                    this.ctx.stroke();
                }
            }
            
            // Disegna la particella con effetto energia turbolenta
            const turbulentSize = particle.size * (1 + particle.turbulentEnergy * 2);
            const turbulentGlow = Math.min(20, particle.turbulentEnergy * 50);
            
            // Glow effect
            if (turbulentGlow > 0) {
                this.ctx.shadowBlur = turbulentGlow;
                this.ctx.shadowColor = `rgba(${particle.color.r}, ${particle.color.g}, ${particle.color.b}, ${alpha * 0.8})`;
            }
            
            this.ctx.fillStyle = `rgba(${particle.color.r}, ${particle.color.g}, ${particle.color.b}, ${alpha})`;
            this.ctx.beginPath();
            this.ctx.arc(particle.x, particle.y, turbulentSize, 0, Math.PI * 2);
            this.ctx.fill();
            
            this.ctx.shadowBlur = 0;
        });
    }
    
    drawTurbulentEnergyField() {
        // Visualizza l'energia cinetica turbolenta come mappa di calore
        const imageData = this.ctx.createImageData(this.width, this.height);
        const data = imageData.data;
        
        for (let y = 0; y < this.height; y += 2) {
            for (let x = 0; x < this.width; x += 2) {
                const col = Math.floor(x / this.cellWidth);
                const row = Math.floor(y / this.cellHeight);
                
                if (col < this.cols && row < this.rows) {
                    const k = this.turbulentKineticEnergy[row][col];
                    const intensity = Math.min(255, k * 300);
                    
                    const index = (y * this.width + x) * 4;
                    data[index] = intensity * 0.3;     // R
                    data[index + 1] = intensity * 0.6; // G
                    data[index + 2] = intensity;       // B
                    data[index + 3] = intensity * 0.15; // A
                }
            }
        }
        
        this.ctx.putImageData(imageData, 0, 0);
    }
    
    
    animate() {
        // Sfondo scuro con leggero fade
        this.ctx.fillStyle = this.colorPalette.background;
        this.ctx.fillRect(0, 0, this.width, this.height);
        
        // Aggiorna la simulazione CFD
        this.updateVortexDynamics();
        this.updateVelocityField();
        this.updateParticles();
        
        // Rendering stratificato per effetto visivo ottimale
        this.drawTurbulentEnergyField();
        this.drawVelocityField();
        this.drawVortices();
        this.drawEddies();
        this.drawParticles();
        
        this.time++;
        requestAnimationFrame(() => this.animate());
    }
    
    resize(width, height) {
        this.width = width;
        this.height = height;
        this.canvas.width = width;
        this.canvas.height = height;
        this.cellWidth = this.width / this.cols;
        this.cellHeight = this.height / this.rows;
        this.initFields();
        this.initVortices();
        this.initEddies();
    }
}

// Inizializzazione migliorata
document.addEventListener('DOMContentLoaded', function() {
    const canvas = document.getElementById('cfd-canvas');
    if (canvas) {
        const container = canvas.parentElement;
        
        function resizeCanvas() {
            const rect = container.getBoundingClientRect();
            canvas.width = rect.width;
            canvas.height = Math.min(500, rect.width * 0.6);
            canvas.style.width = '100%';
            canvas.style.height = canvas.height + 'px';
        }
        
        resizeCanvas();
        const cfd = new AdvancedCFDAnimation(canvas);
        
        window.addEventListener('resize', () => {
            resizeCanvas();
            cfd.resize(canvas.width, canvas.height);
        });
        
        // Pausa l'animazione quando la pagina non è visibile (performance)
        document.addEventListener('visibilitychange', function() {
            if (document.hidden) {
                // L'animazione continuerà ma con meno dettagli se necessario
            }
        });
    }
});
