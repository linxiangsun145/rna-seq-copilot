import type { Metadata } from "next";
import "./globals.css";
import { Toaster } from "sonner";

export const metadata: Metadata = {
  title: "RNA-seq Copilot",
  description: "Intelligent RNA-seq count matrix analysis and differential expression assistant",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" translate="no" className="notranslate" suppressHydrationWarning>
      <head>
        <meta name="google" content="notranslate" />
        <meta httpEquiv="Content-Language" content="en" />
      </head>
      <body className="antialiased" suppressHydrationWarning>
        <div className="min-h-screen bg-gradient-to-br from-slate-50 to-blue-50">
          <nav className="border-b bg-white/80 backdrop-blur-sm sticky top-0 z-50" suppressHydrationWarning>
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
              <div className="flex items-center justify-between h-16">
                <a href="/" className="flex items-center gap-2">
                  <div className="w-8 h-8 rounded-lg bg-blue-600 flex items-center justify-center">
                    <span className="text-white font-bold text-sm" suppressHydrationWarning>R</span>
                  </div>
                  <span className="font-semibold text-gray-900" suppressHydrationWarning>RNA-seq Copilot</span>
                </a>
                <div className="flex items-center gap-6 text-sm text-gray-600">
                  <a href="/" className="hover:text-blue-600 transition-colors">Home</a>
                  <a href="/upload" className="hover:text-blue-600 transition-colors">Upload</a>
                  <a href="/about" className="hover:text-blue-600 transition-colors">About</a>
                  <a href="/help" className="hover:text-blue-600 transition-colors">Help</a>
                </div>
              </div>
            </div>
          </nav>
          <main suppressHydrationWarning>{children}</main>
        </div>
        <Toaster position="top-right" richColors />
      </body>
    </html>
  );
}
