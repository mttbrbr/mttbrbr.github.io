import { getCollection, type CollectionEntry } from 'astro:content';
import katex from 'katex';

export type Kind = 'project' | 'note';

/** Personal notes are hidden for now: flip to true to publish them again. */
export const SHOW_NOTES = false;

export interface Entry {
  kind: Kind;
  id: string;
  url: string;
  title: string;
  date: Date;
  tags: string[];
  excerpt: string;
  excerptHtml: string;
  minutes: number;
  searchText: string;
}

/** Render inline $...$ math in short strings (excerpts) at build time. */
export function inlineMath(text: string): string {
  const escape = (s: string) =>
    s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  return text
    .split(/(\$[^$]+\$)/g)
    .map((part) =>
      part.startsWith('$') && part.endsWith('$') && part.length > 2
        ? katex.renderToString(part.slice(1, -1), { throwOnError: false })
        : escape(part),
    )
    .join('');
}

function plain(md: string): string {
  return md
    .replace(/\$\$[\s\S]*?\$\$/g, ' ')
    .replace(/```[\s\S]*?```/g, ' ')
    .replace(/<[^>]+>/g, ' ')
    .replace(/[#*_`>\[\]()|!-]/g, ' ')
    .replace(/\s+/g, ' ')
    .toLowerCase();
}

function build(
  kind: Kind,
  e: CollectionEntry<'projects'> | CollectionEntry<'notes'>,
  tags: string[],
): Entry {
  const body = e.body ?? '';
  const words = body.split(/\s+/).length;
  return {
    kind,
    id: e.id,
    url: `/${kind}s/${e.id}/`,
    title: e.data.title,
    date: e.data.date,
    tags,
    excerpt: e.data.excerpt,
    excerptHtml: inlineMath(e.data.excerpt),
    minutes: Math.max(1, Math.round(words / 220)),
    searchText: [e.data.title, e.data.excerpt, tags.join(' '), plain(body)]
      .join(' ')
      .toLowerCase(),
  };
}

/** Projects and notes merged into a single, date-sorted feed. */
export async function getAllEntries(): Promise<Entry[]> {
  const projects = await getCollection('projects');
  const notes = SHOW_NOTES ? await getCollection('notes', (n) => n.data.published) : [];
  return [
    ...projects.map((p) => build('project', p, p.data.tech)),
    ...notes.map((n) => build('note', n, n.data.categories ? [n.data.categories] : [])),
  ].sort((a, b) => b.date.valueOf() - a.date.valueOf());
}

export const fmtDate = (d: Date) =>
  d.toLocaleDateString('en-GB', { year: 'numeric', month: 'short', day: '2-digit' });
