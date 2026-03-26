"use client";

async function downloadBlob(url: string, filename: string) {
  try {
    const res = await fetch(url);
    const blob = await res.blob();
    const blobUrl = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = blobUrl;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(blobUrl);
  } catch {
    window.open(url, "_blank");
  }
}

interface Props {
  src: string;
  title: string;
  description?: string;
}

export function PlotViewer({ src, title, description }: Props) {
  return (
    <div className="rounded-xl border bg-white overflow-hidden">
      <div className="px-4 py-3 border-b">
        <h3 className="font-medium text-gray-900 text-sm">{title}</h3>
        {description && (
          <p className="text-xs text-gray-500 mt-0.5">{description}</p>
        )}
      </div>
      <div className="p-4 flex items-center justify-center bg-gray-50 min-h-[280px]">
        {/* eslint-disable-next-line @next/next/no-img-element */}
        <img
          src={src}
          alt={title}
          className="max-w-full max-h-[400px] object-contain rounded"
          loading="lazy"
        />
      </div>
      <div className="px-4 py-2 border-t bg-white">
        <button
          onClick={() => downloadBlob(src, src.split("/").pop() || "plot.png")}
          className="text-xs text-blue-600 hover:underline"
        >
          Download PNG
        </button>
      </div>
    </div>
  );
}
