"use client";

import type { DEGRow } from "@/lib/types";
import { useState } from "react";

interface Props {
  rows: DEGRow[];
  downloadUrl: string;
}

const PAGE_SIZE = 20;

export function DEGTable({ rows, downloadUrl }: Props) {
  const [page, setPage] = useState(0);
  const [sortKey, setSortKey] = useState<keyof DEGRow>("padj");
  const [sortAsc, setSortAsc] = useState(true);
  const [filter, setFilter] = useState("");

  const filtered = rows.filter((r) =>
    r.gene_id.toLowerCase().includes(filter.toLowerCase())
  );

  const sorted = [...filtered].sort((a, b) => {
    const va = a[sortKey];
    const vb = b[sortKey];
    if (typeof va === "string" && typeof vb === "string") {
      return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    }
    const na = (va as number) ?? NaN;
    const nb = (vb as number) ?? NaN;
    if (isNaN(na) && isNaN(nb)) return 0;
    if (isNaN(na)) return 1;
    if (isNaN(nb)) return -1;
    return sortAsc ? na - nb : nb - na;
  });

  const totalPages = Math.ceil(sorted.length / PAGE_SIZE);
  const paged = sorted.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE);

  const handleSort = (key: keyof DEGRow) => {
    if (key === sortKey) setSortAsc((v) => !v);
    else { setSortKey(key); setSortAsc(true); }
    setPage(0);
  };

  const col = (key: keyof DEGRow, label: string) => (
    <th
      key={key}
      className="px-3 py-2 text-left text-xs font-medium text-gray-500 uppercase cursor-pointer hover:text-gray-900 select-none"
      onClick={() => handleSort(key)}
    >
      {label} {sortKey === key ? (sortAsc ? "↑" : "↓") : ""}
    </th>
  );

  const fmt = (v: number, digits = 3) =>
    v === null || v === undefined || isNaN(v) ? "NA" : v.toFixed(digits);

  const fmtSci = (v: number) =>
    v === null || isNaN(v) ? "NA" : v.toExponential(2);

  return (
    <div className="rounded-xl border bg-white overflow-hidden">
      <div className="flex items-center justify-between px-4 py-3 border-b gap-3">
        <h3 className="font-medium text-gray-900 text-sm">
          DEG Results ({filtered.length.toLocaleString()} genes)
        </h3>
        <div className="flex items-center gap-2">
          <input
            type="text"
            placeholder="Search gene…"
            value={filter}
            onChange={(e) => { setFilter(e.target.value); setPage(0); }}
            className="border rounded px-2 py-1 text-xs w-36 focus:outline-none focus:ring-1 focus:ring-blue-400"
          />
          <a
            href={downloadUrl}
            download
            className="px-3 py-1 bg-blue-600 text-white text-xs rounded hover:bg-blue-700 transition-colors"
          >
            Download CSV
          </a>
        </div>
      </div>

      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead className="bg-gray-50 border-b">
            <tr>
              {col("gene_id", "Gene")}
              {col("baseMean", "baseMean")}
              {col("log2FoldChange", "log2FC")}
              {col("pvalue", "p-value")}
              {col("padj", "padj")}
            </tr>
          </thead>
          <tbody className="divide-y divide-gray-100">
            {paged.map((row) => (
              <tr key={row.gene_id} className="hover:bg-gray-50">
                <td className="px-3 py-2 font-mono text-xs text-gray-900">{row.gene_id}</td>
                <td className="px-3 py-2 text-gray-600 text-xs">{fmt(row.baseMean, 1)}</td>
                <td
                  className={`px-3 py-2 text-xs font-medium ${
                    row.log2FoldChange > 0 ? "text-red-600" : "text-blue-600"
                  }`}
                >
                  {fmt(row.log2FoldChange)}
                </td>
                <td className="px-3 py-2 text-gray-600 text-xs font-mono">{fmtSci(row.pvalue)}</td>
                <td className="px-3 py-2 text-xs font-mono">
                  <span
                    className={
                      row.padj < 0.05 ? "text-red-600 font-semibold" : "text-gray-600"
                    }
                  >
                    {row.padj ? fmtSci(row.padj) : "NA"}
                  </span>
                </td>
              </tr>
            ))}
            {paged.length === 0 && (
              <tr>
                <td colSpan={5} className="px-3 py-6 text-center text-xs text-gray-400">
                  No genes match the filter.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {totalPages > 1 && (
        <div className="flex items-center justify-between px-4 py-2 border-t text-xs text-gray-500">
          <span>
            Page {page + 1} of {totalPages}
          </span>
          <div className="flex gap-2">
            <button
              onClick={() => setPage((p) => Math.max(0, p - 1))}
              disabled={page === 0}
              className="px-2 py-1 rounded border disabled:opacity-40 hover:bg-gray-50"
            >
              Prev
            </button>
            <button
              onClick={() => setPage((p) => Math.min(totalPages - 1, p + 1))}
              disabled={page === totalPages - 1}
              className="px-2 py-1 rounded border disabled:opacity-40 hover:bg-gray-50"
            >
              Next
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
