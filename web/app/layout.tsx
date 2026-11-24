import React from 'react'
import localFont from 'next/font/local'
import './globals.css'

// Self-hosted, OFL 1.1 (licences next to the files in app/fonts).
const text = localFont({ src: './fonts/AtkinsonHyperlegibleNext-wght.woff2', weight: '200 800', variable: '--font-text', display: 'swap' })
const mono = localFont({ src: './fonts/JetBrainsMono-wght.woff2', weight: '100 800', variable: '--font-mono', display: 'swap' })

export const metadata = {
  title: 'Parallel Pangenome Graph Explorer',
  description: 'Comparing monolithic and parallel/reassembled pangenome graph construction. BCM SV Hackathon 2026.',
  openGraph: {
    title: 'Parallel Pangenome Graph Explorer',
    description: 'Can regional pangenome graphs be built in parallel and stitched back together?',
    type: 'website',
  },
}

// Matches --paper in globals.css; meta tags need a literal colour.
export const viewport = { themeColor: '#f8fafb', colorScheme: 'light' }

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en" className={`${text.variable} ${mono.variable}`}>
      <body>{children}</body>
    </html>
  )
}
