---
layout: page
title: Contatti
permalink: /contatti/
---

# Contatti

Sono sempre interessato a discutere nuove opportunità, progetti e collaborazioni.

<div class="contact-grid">
    <div class="contact-info card">
        <h2>Info di contatto</h2>
        
        <ul class="contact-list">
            <li>
                <strong>Email:</strong> 
                <a href="mailto:{{ site.email }}">{{ site.email }}</a>
            </li>
            <li>
                <strong>LinkedIn:</strong> 
                <a href="https://linkedin.com/in/marcorossi-dev" target="_blank">linkedin.com/in/marcorossi-dev</a>
            </li>
            <li>
                <strong>GitHub:</strong> 
                <a href="https://github.com/{{ site.github_username }}" target="_blank">github.com/{{ site.github_username }}</a>
            </li>
        </ul>
    </div>
    
    <div class="contact-form card">
        <h2>Invia un messaggio</h2>
        <p>Compila il form sottostante per contattarmi direttamente.</p>
        
        <form id="contact-form" action="https://formspree.io/f/your-formspree-id" method="POST">
            <div class="form-group">
                <label for="name">Nome</label>
                <input type="text" id="name" name="name" required>
            </div>
            
            <div class="form-group">
                <label for="email">Email</label>
                <input type="email" id="email" name="email" required>
            </div>
            
            <div class="form-group">
                <label for="subject">Oggetto</label>
                <input type="text" id="subject" name="subject" required>
            </div>
            
            <div class="form-group">
                <label for="message">Messaggio</label>
                <textarea id="message" name="message" rows="5" required></textarea>
            </div>
            
            <button type="submit" class="btn btn-primary">Invia messaggio</button>
        </form>
    </div>
</div>

## Orari di disponibilità

Sono generalmente disponibile per call e meeting dal lunedì al venerdì, dalle 9:00 alle 18:00 (CET).

Per richieste urgenti, puoi contattarmi direttamente via email indicando [URGENTE] nell'oggetto.
