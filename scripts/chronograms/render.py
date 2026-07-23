#!/usr/bin/env python3

import sys
import re
import html
import argparse
import gzip
import bz2
import lzma
import datetime
import json
import colorsys

try:
    import compression.zstd
    has_zstd = True
except ModuleNotFoundError:
    has_zstd = False


def feed(files):
    if not files:
        for c in sys.stdin.readlines():
            yield c
    else:
        for file in files:
            if re.match(r"^.*\.zstd?$", file):
                assert has_zstd
                with compression.zstd.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*\.gz$", file):
                with gzip.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*\.xz$", file):
                with lzma.open(file, 'r') as f:
                    for c in f.readlines():
                        yield c.decode()
            elif re.match(r"^.*\.bz2$", file):
                with bz2.open(file) as f:
                    for c in f.readlines():
                        yield c.decode()
            else:
                with open(file) as f:
                    for c in f.readlines():
                        yield c


def convert_epoch_to_readable(t):
    dt = datetime.datetime.fromtimestamp(int(t)).astimezone()
    nsec = int((t - int(t)) * 1e9)
    iso_fmt = f"{dt.strftime('%Y-%m-%d %H:%M:%S')}"
    iso_fmt += f".{nsec:09d}{dt.strftime('%z')}"
    return iso_fmt


PREDEFINED_COLORS = {
    "INIT": "rgb(78, 121, 167)",
    "SKEWGAUSS": "rgb(242, 142, 43)",
    "ADJUST": "rgb(225, 87, 89)",
    "ALLOC": "rgb(118, 183, 178)",
    "FIB": "rgb(89, 161, 79)",
    "SSS": "rgb(237, 201, 72)",
    "PBR": "rgb(176, 122, 161)",
    "ECM": "rgb(255, 157, 167)",
    "DS": "rgb(156, 117, 95)",
    "PCLAT": "rgb(186, 176, 172)",
    "SLICING": "rgb(92, 192, 222)",
    "BOTCHED": "rgb(200, 0, 0)",
}


def get_color_str(category):
    cat_upper = category.upper()
    if cat_upper in PREDEFINED_COLORS:
        return PREDEFINED_COLORS[cat_upper]

    hue = (hash(cat_upper) & 0xFFFFFFFF) * 0.618033988749895 % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.85)
    return f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"


LINE_REGEX = re.compile(r"^t\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+(.+)$")


def parse_trace_data(input_stream):
    bubbles = []
    for line in input_stream:
        line = line.strip()
        if not line or not line.startswith('t '):
            continue
        match = LINE_REGEX.match(line)
        if not match:
            continue

        thread_id = int(match.group(1))
        start_time = float(match.group(2))
        duration = float(match.group(3))
        metric = int(match.group(4))
        description = match.group(5).strip()

        category = description.split()[0] if description else "UNKNOWN"

        bubbles.append({
            "thread": thread_id,
            "start": start_time,
            "duration": duration,
            "end": start_time + duration,
            "metric": metric,
            "desc": description,
            "category": category
        })
    return bubbles


def clamp_trace_data(bubbles, tstart_arg, tend_arg):
    if not bubbles:
        return bubbles, 0.0, 0.0

    raw_min_start = min(b["start"] for b in bubbles)
    raw_max_end = max(b["end"] for b in bubbles)

    if tstart_arg is None:
        win_start = raw_min_start
    elif tstart_arg > 1e8:
        win_start = tstart_arg
    else:
        win_start = raw_min_start + tstart_arg

    if tend_arg is None:
        win_end = raw_max_end
    elif tend_arg > 1e8:
        win_end = tend_arg
    else:
        win_end = raw_min_start + tend_arg

    if win_start >= win_end:
        print(f"Warning: Clamp range [{win_start}, {win_end}] is invalid.", file=sys.stderr)
        return [], win_start, win_end

    clamped = []
    for b in bubbles:
        if b["end"] <= win_start or b["start"] >= win_end:
            continue

        cb = dict(b)
        cb["start"] = max(b["start"], win_start)
        cb["end"] = min(b["end"], win_end)
        cb["duration"] = cb["end"] - cb["start"]
        clamped.append(cb)

    return clamped, win_start, win_end


def generate_xhtml_chronogram(bubbles, win_start, win_end, output_filepath="chronogram.html"):
    if not bubbles:
        print("Error: No valid time bubble data found in window.", file=sys.stderr)
        return

    threads = sorted(list({b["thread"] for b in bubbles}))
    thread_map = {t: idx for idx, t in enumerate(threads)}

    # Build unique categories list and assign color indices
    categories = sorted(list({b["category"] for b in bubbles}))
    category_map = {cat: idx for idx, cat in enumerate(categories)}
    category_colors = [get_color_str(cat) for cat in categories]

    # Compact array format: { col_idx: [ [start, duration, cat_idx, metric, desc], ... ] }
    trace_data_js = {}
    for b in bubbles:
        col_idx = thread_map[b["thread"]]
        if col_idx not in trace_data_js:
            trace_data_js[col_idx] = []

        trace_data_js[col_idx].append([
            round(b["start"] - win_start, 9),
            round(b["duration"], 9),
            category_map[b["category"]],
            b["metric"],
            b["desc"]
        ])

    # Ensure events in each column are sorted by start time
    for col_idx in trace_data_js:
        trace_data_js[col_idx].sort(key=lambda x: x[0])

    json_trace_data = json.dumps(trace_data_js, separators=(',', ':'))

    legend_html = []
    for cat in categories:
        color = get_color_str(cat)
        legend_html.append(
            f'<span class="legend-item"><span class="color-box" style="background-color: {color};"></span>{cat}</span>'
        )

    xhtml_content = f"""<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>cado-nfs chronograms</title>
<style type="text/css">
  html, body {{
    height: 100%;
    margin: 0;
    padding: 0;
    overflow: hidden;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
    background-color: #f0f0f0;
  }}

  .app-container {{
    display: flex;
    flex-direction: column;
    height: 100vh;
    padding: 12px;
    box-sizing: border-box;
  }}

  .top-frame {{
    flex: 0 0 auto;
    background: #ffffff;
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 12px 16px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.06);
    margin-bottom: 10px;
  }}

  .header-row {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 6px;
  }}

  .header-title {{
    margin: 0;
    font-size: 16px;
    color: #222;
  }}

  .controls {{
    font-size: 12px;
    color: #444;
    display: flex;
    align-items: center;
    gap: 12px;
    background: #f0f0f0;
    padding: 4px 10px;
    border-radius: 4px;
    border: 1px solid #e0e0e0;
  }}

  .controls label {{
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 3px;
  }}

  .window-meta {{
    font-size: 11px;
    color: #666;
    margin-bottom: 8px;
  }}

  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin-bottom: 8px;
    padding-bottom: 8px;
    border-bottom: 1px solid #eee;
  }}

  .legend-item {{
    display: inline-flex;
    align-items: center;
    font-size: 12px;
    color: #444;
  }}

  .color-box {{
    width: 12px;
    height: 12px;
    display: inline-block;
    margin-right: 5px;
    border: 1px solid rgba(0,0,0,0.15);
  }}

  #information {{
    margin: 0;
    padding: 8px 12px;
    background: #f8f9fa;
    color: #333333;
    border: 1px solid #dcdcdc;
    border-radius: 3px;
    font-family: SFMono-Regular, Consolas, "Liberation Mono", Menlo, monospace;
    font-size: 12px;
    min-height: 18px;
    white-space: pre-wrap;
    word-break: break-all;
  }}

  .chart-viewport {{
    flex: 1 1 auto;
    position: relative;
    background: #ffffff;
    border: 1px solid #bbb;
    border-radius: 4px;
    overflow: hidden;
  }}

  #viewport-canvas {{
    width: 100%;
    height: 100%;
    display: block;
    cursor: crosshair;
  }}
</style>
</head>

<body>

<div class="app-container">

  <div class="top-frame">
    <div class="header-row">
      <h3 class="header-title">cado-nfs Execution Chronogram</h3>

      <div class="controls">
        <strong>Delimiters:</strong>
        <label><input type="radio" name="border-opt" value="none" checked="checked" /> Off</label>
        <label><input type="radio" name="border-opt" value="top" /> On</label>
        <span style="color:#aaa;">|</span>
        <span style="font-size: 11px; color: #666;">
          <strong>Scroll:</strong> Pan Time &nbsp;•&nbsp; <strong>Ctrl + Scroll:</strong> Zoom Time
        </span>
      </div>
    </div>

    <div class="window-meta">
      Time Window:
      <strong>{convert_epoch_to_readable(win_start)}</strong> to
      <strong>{convert_epoch_to_readable(win_end)}</strong>
      (Span: {(win_end-win_start)*1000:.3f} ms)
      | Rendered Events: {len(bubbles)}
    </div>

    <div class="legend">
      {' '.join(legend_html)}
    </div>

    <div id="information">Initializing viewport canvas...</div>
  </div>

  <div class="chart-viewport">
    <canvas id="viewport-canvas"></canvas>
  </div>

</div>

<script type="text/javascript">
//<![CDATA[
const TRACE = {json_trace_data};
const THREADS = {threads};
const CATEGORY_NAMES = {json.dumps(categories)};
const CATEGORY_COLORS = {json.dumps(category_colors)};

const TOTAL_SPAN = {win_end - win_start};

document.addEventListener("DOMContentLoaded", () => {{
  const infoElement = document.getElementById("information");
  const canvas = document.getElementById("viewport-canvas");
  const ctx = canvas.getContext("2d");

  const LEFT_MARGIN = 75;
  const RIGHT_MARGIN = 20;
  const HEADER_HEIGHT = 30;
  const BOTTOM_MARGIN = 20;

  let showDelimiters = false;
  let activeBubble = null;

  // Viewport State (in relative seconds)
  let viewStart = 0.0;
  let viewSpan = TOTAL_SPAN;

  function resizeCanvas() {{
    const rect = canvas.parentElement.getBoundingClientRect();
    canvas.width = rect.width;
    canvas.height = rect.height;
    render();
  }}

  // Binary search to find index of first event ending after targetStart
  function findFirstVisibleIndex(events, targetStart) {{
    let low = 0, high = events.length - 1, ans = events.length;
    while (low <= high) {{
      const mid = (low + high) >> 1;
      if (events[mid][0] + events[mid][1] >= targetStart) {{
        ans = mid;
        high = mid - 1;
      }} else {{
        low = mid + 1;
      }}
    }}
    return ans;
  }}

  function render() {{
    const W = canvas.width;
    const H = canvas.height;
    ctx.clearRect(0, 0, W, H);

    const numThreads = THREADS.length;
    const perThreadWidth = Math.max(4, Math.floor((W - LEFT_MARGIN - RIGHT_MARGIN) / numThreads));
    const renderWidth = numThreads * perThreadWidth;

    const availableH = H - HEADER_HEIGHT - BOTTOM_MARGIN;
    const y_scale = availableH / viewSpan;
    const viewEnd = viewStart + viewSpan;

    // 1. Draw Time Grid Lines & Labels
    const numTicks = 8;
    ctx.fillStyle = "#777777";
    ctx.strokeStyle = "#ebebeb";
    ctx.font = "11px sans-serif";
    ctx.lineWidth = 1;

    for (let i = 0; i <= numTicks; i++) {{
      const rel_t = viewStart + (viewSpan * (i / numTicks));
      const y = HEADER_HEIGHT + Math.floor(i * (availableH / numTicks));

      ctx.beginPath();
      ctx.moveTo(LEFT_MARGIN - 5, y);
      ctx.lineTo(LEFT_MARGIN + renderWidth, y);
      ctx.stroke();

      const ms = rel_t * 1000;
      const msStr = Math.abs(ms) < 0.005 ? "0.00" : ms.toFixed(2);
      const label = (msStr.startsWith('-') ? "" : "+") + msStr;
      ctx.fillText(`${label}ms`, 5, y + 4);
    }}

    // 2. Draw Column Headers & Vertical Dividers
    ctx.fillStyle = "#222222";
    ctx.strokeStyle = "#d0d0d0";
    for (let i = 0; i < numThreads; i++) {{
      const x1 = LEFT_MARGIN + i * perThreadWidth;
      ctx.fillText(`T${{THREADS[i]}}`, x1 + 4, 18);

      ctx.beginPath();
      ctx.moveTo(x1, HEADER_HEIGHT);
      ctx.lineTo(x1, H - BOTTOM_MARGIN);
      ctx.stroke();
    }}

    // 3. Render Visible Bubbles using Binary Search
    for (let colIdx = 0; colIdx < numThreads; colIdx++) {{
      const events = TRACE[colIdx];
      if (!events) continue;

      const x1 = LEFT_MARGIN + colIdx * perThreadWidth;
      const startIdx = findFirstVisibleIndex(events, viewStart);

      for (let i = startIdx; i < events.length; i++) {{
        const e = events[i];
        const eStart = e[0];
        const eEnd = e[0] + e[1];

        if (eStart > viewEnd) break; // Past the visible bottom

        const y1 = HEADER_HEIGHT + (eStart - viewStart) * y_scale;
        const y2 = HEADER_HEIGHT + (eEnd - viewStart) * y_scale;
        const h = Math.max(1, y2 - y1);

        // Draw Bubble
        ctx.fillStyle = CATEGORY_COLORS[e[2]];
        ctx.fillRect(x1, y1, perThreadWidth, h);

        // Optional Top Delimiter
        if (showDelimiters) {{
          ctx.fillStyle = "rgba(0, 0, 0, 0.6)";
          ctx.fillRect(x1, y1, perThreadWidth, 1);
        }}
      }}
    }}

    // 4. Draw Active Hover Highlight
    if (activeBubble) {{
      const [colIdx, e] = activeBubble;
      const x1 = LEFT_MARGIN + colIdx * perThreadWidth;
      const y1 = HEADER_HEIGHT + (e[0] - viewStart) * y_scale;
      const y2 = HEADER_HEIGHT + ((e[0] + e[1]) - viewStart) * y_scale;
      const h = Math.max(1, y2 - y1);

      ctx.fillStyle = "rgba(0, 0, 0, 0.35)";
      ctx.fillRect(x1, y1, perThreadWidth, h);
      ctx.strokeStyle = "#000000";
      ctx.strokeRect(x1, y1, perThreadWidth, h);
    }}
  }}

  // --- Interaction: Hover Detection ---
  canvas.addEventListener("mousemove", (e) => {{
    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    const numThreads = THREADS.length;
    const perThreadWidth = Math.max(4, Math.floor((canvas.width - LEFT_MARGIN - RIGHT_MARGIN) / numThreads));

    const colIdx = Math.floor((mouseX - LEFT_MARGIN) / perThreadWidth);
    let hit = null;

    if (colIdx >= 0 && colIdx < numThreads && mouseY >= HEADER_HEIGHT && mouseY <= canvas.height - BOTTOM_MARGIN) {{
      const availableH = canvas.height - HEADER_HEIGHT - BOTTOM_MARGIN;
      const tMouse = viewStart + ((mouseY - HEADER_HEIGHT) / availableH) * viewSpan;

      const events = TRACE[colIdx];
      if (events) {{
        const idx = findFirstVisibleIndex(events, tMouse);
        if (idx < events.length && events[idx][0] <= tMouse) {{
          hit = [colIdx, events[idx]];
        }}
      }}
    }}

    if (hit !== activeBubble) {{
      activeBubble = hit;
      if (hit) {{
        const e = hit[1];
        const catName = CATEGORY_NAMES[e[2]];
        infoElement.textContent = `Thread ${{THREADS[hit[0]]}} | ${{catName}} | Start: +${{(e[0]*1000).toFixed(3)}}ms | Duration: ${{(e[1]*1000).toFixed(3)}}ms | Metric: ${{e[3]}} | Details: ${{e[4]}}`;
      }} else {{
        infoElement.textContent = "Hover over a region to inspect details...";
      }}
      render();
    }}
  }});

  canvas.addEventListener("mouseleave", () => {{
    activeBubble = null;
    infoElement.textContent = "Hover over a region to inspect details...";
    render();
  }});

  // --- Interaction: Scroll (Pan) and Ctrl+Scroll (Zoom) ---
  canvas.addEventListener("wheel", (e) => {{
    e.preventDefault();

    if (e.ctrlKey || e.metaKey) {{
      // ZOOM at mouse Y position
      const rect = canvas.getBoundingClientRect();
      const mouseY = e.clientY - rect.top;
      const availableH = canvas.height - HEADER_HEIGHT - BOTTOM_MARGIN;
      const mouseRelT = Math.max(0, Math.min(1, (mouseY - HEADER_HEIGHT) / availableH));
      const tFocus = viewStart + mouseRelT * viewSpan;

      const zoomFactor = e.deltaY > 0 ? 1.15 : 0.85;
      const newSpan = Math.max(0.000001, Math.min(TOTAL_SPAN * 2, viewSpan * zoomFactor));

      viewStart = tFocus - mouseRelT * newSpan;
      viewSpan = newSpan;
    }} else {{
      // PAN vertically
      const panAmount = (e.deltaY / canvas.height) * viewSpan * 0.5;
      viewStart = Math.max(0.0 - viewSpan * 0.1, Math.min(TOTAL_SPAN - viewSpan * 0.9, viewStart + panAmount));
    }}

    render();
  }}, {{ passive: false }});

  // Border option listeners
  document.querySelectorAll("input[name='border-opt']").forEach(radio => {{
    radio.addEventListener("change", (e) => {{
      showDelimiters = (e.target.value === "top");
      render();
    }});
  }});

  window.addEventListener("resize", resizeCanvas);
  resizeCanvas();
  infoElement.textContent = "Hover over a region to inspect details...";
}});
//]]>
</script>

</body>
</html>
"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(xhtml_content)

    print(f"Generated Canvas Chronogram: {output_filepath}")
    print(f"Time Range: {convert_epoch_to_readable(win_start)} -> {convert_epoch_to_readable(win_end)}")


def main():
    parser = argparse.ArgumentParser(description="Generate an interactive Canvas Chronogram from CADO-NFS trace logs.")
    parser.add_argument("file", nargs='*', help="Path to trace file (default: stdin)")
    parser.add_argument("-o", "--output", default="chronogram.html", help="Output HTML path (default: chronogram.html)")
    parser.add_argument("--tstart", type=float, default=None, help="Start time offset or absolute timestamp.")
    parser.add_argument("--tend", type=float, default=None, help="End time offset or absolute timestamp.")

    args = parser.parse_args()

    raw_data = parse_trace_data(feed(args.file))
    clamped_data, win_start, win_end = clamp_trace_data(raw_data, args.tstart, args.tend)

    generate_xhtml_chronogram(clamped_data, win_start, win_end, args.output)


if __name__ == "__main__":
    main()
