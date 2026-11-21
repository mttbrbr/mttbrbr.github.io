# mttbrbr.com

Personal site built with [Astro](https://astro.build), deployed to GitHub Pages.

```sh
npm install
npm run dev       # http://localhost:4321, live reload
npm run build     # static output in dist/
```

## Adding content

- **Project:** add `src/content/projects/<slug>.md` (or `.mdx` for interactive components)
  with `title`, `date`, `tech`, `excerpt` in the frontmatter.
- **Note:** add `src/content/notes/<slug>.md` with `title`, `date`, `categories`, `excerpt`.
- Images and files go in `public/` and are referenced as `/path/to/file.png`.

Both collections are merged into the searchable index on the homepage. LaTeX (`$...$`, `$$...$$`) is rendered at build time with KaTeX.
