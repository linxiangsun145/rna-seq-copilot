"use client";

import { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import { cn } from "@/lib/utils";

interface FileUploaderProps {
  label: string;
  accept?: string;
  hint?: string;
  onFile: (file: File) => void;
  file: File | null;
}

export function FileUploader({
  label,
  accept = ".csv,.tsv,.txt",
  hint,
  onFile,
  file,
}: FileUploaderProps) {
  const [hover, setHover] = useState(false);

  const onDrop = useCallback(
    (accepted: File[]) => {
      if (accepted[0]) onFile(accepted[0]);
    },
    [onFile]
  );

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: { "text/plain": [".csv", ".tsv", ".txt"], "text/csv": [".csv"] },
    maxFiles: 1,
  });

  return (
    <div className="space-y-1">
      <label className="text-sm font-medium text-gray-700">{label}</label>
      <div
        {...getRootProps()}
        className={cn(
          "border-2 border-dashed rounded-lg p-6 text-center cursor-pointer transition-colors",
          isDragActive
            ? "border-blue-400 bg-blue-50"
            : file
            ? "border-green-400 bg-green-50"
            : "border-gray-200 hover:border-blue-300 hover:bg-blue-50/40"
        )}
      >
        <input {...getInputProps()} />
        {file ? (
          <div className="flex items-center justify-center gap-2">
            <span className="text-green-600 text-lg">✓</span>
            <span className="text-sm font-medium text-green-700">{file.name}</span>
            <span className="text-xs text-gray-400">
              ({(file.size / 1024).toFixed(1)} KB)
            </span>
          </div>
        ) : (
          <div>
            <div className="text-gray-400 text-2xl mb-2">📄</div>
            <p className="text-sm text-gray-600">
              {isDragActive ? "Drop here…" : "Drag & drop or click to browse"}
            </p>
            {hint && <p className="text-xs text-gray-400 mt-1">{hint}</p>}
          </div>
        )}
      </div>
    </div>
  );
}
