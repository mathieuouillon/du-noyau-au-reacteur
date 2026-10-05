---
title: Bibliographie
icon: fas fa-book-open
order: 3
---

{% include cours-lecons.html %}
{% assign biblio = site.data.bibliographie %}

Toutes les références du cours, regroupées par type. Chaque leçon donne, en fin
de page, sa propre bibliographie commentée et la source de ses chiffres. Les
DOI ont été vérifiés un par un auprès de [doi.org](https://doi.org) et de
[Crossref](https://www.crossref.org).

**Pour commencer :** en français, P. Reuss, *Précis de neutronique*, et
L. Valentin, *Physique subatomique* ; en anglais, K.S. Krane, *Introductory
Nuclear Physics*. Pour la fission, R. Vandenbosch et J.R. Huizenga, *Nuclear
Fission*, puis H.J. Krappe et K. Pomorski, *Theory of Nuclear Fission*.

{% assign groupes = "livre|Manuels et ouvrages,revue|Articles de synthèse,article|Articles,donnees|Données évaluées,code|Codes et outils" | split: "," %}
{% for g in groupes %}
{% assign gt = g | split: "|" %}
<h2 id="{{ gt[0] }}">{{ gt[1] }}</h2>
<ul class="biblio">
{%- for paire in biblio -%}
  {%- assign cle = paire[0] -%}
  {%- assign e = paire[1] -%}
  {%- if e.type == gt[0] %}
  <li id="ref-{{ cle }}">{% include biblio-entree.html e=e %}
    {%- assign cite = "" -%}
    {%- for l in cours_lecons -%}
      {%- for item in l.bibliographie -%}
        {%- assign c = item.cle | default: item -%}
        {%- if c == cle -%}
          {%- capture lien -%}<a href="{{ l.url | relative_url }}">{% if l.annexe %}annexe{% else %}{{ l.lecon }}{% endif %}</a>{%- endcapture -%}
          {%- assign cite = cite | append: lien | append: " " -%}
        {%- endif -%}
      {%- endfor -%}
    {%- endfor -%}
    {%- if cite != "" %}<br><span class="biblio-note">Cité dans la leçon : {{ cite | strip | replace: "> <", ">, <" }}</span>{% endif %}
  </li>
  {%- endif -%}
{%- endfor %}
</ul>
{% endfor %}
