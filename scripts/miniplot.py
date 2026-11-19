#!/usr/bin/env python3
"""Tiny dependency-free plotting helpers (Pillow only) used for the figures.

This is deliberately small: axes, ticks, lines, scatter, images, legends and
vertical stacking.  It keeps the repository reproducible without matplotlib.
"""
from PIL import Image, ImageDraw, ImageFont

FONT_PATHS = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
]


def font(size, bold=False):
    for path in FONT_PATHS:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _nice_ticks(lo, hi, target=6):
    if hi <= lo:
        return [lo]
    span = hi - lo
    raw = span / max(1, target - 1)
    import math
    mag = 10 ** math.floor(math.log10(raw))
    norm = raw / mag
    if norm < 1.5:
        step = 1.0
    elif norm < 3:
        step = 2.0
    elif norm < 7:
        step = 5.0
    else:
        step = 10.0
    step *= mag
    start = math.ceil(lo / step) * step
    ticks = []
    v = start
    while v <= hi + step * 1e-6:
        ticks.append(0.0 if abs(v) < step * 1e-9 else v)
        v += step
    return ticks


def _fmt(v):
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1e4 or av < 1e-3:
        return f"{v:.1e}"
    if abs(v - round(v)) < 1e-9:
        return str(int(round(v)))
    return f"{v:.3g}"


class Canvas:
    def __init__(self, width=1400, height=820, scale=2,
                 margins=(95, 45, 55, 85), background="white"):
        self.w = width
        self.h = height
        self.s = scale
        self.W = width * scale
        self.H = height * scale
        self.img = Image.new("RGB", (self.W, self.H), background)
        self.draw = ImageDraw.Draw(self.img)
        left, right, top, bottom = margins
        self.rect = (left * scale, top * scale,
                     (width - right) * scale, (height - bottom) * scale)
        self.xlim = None
        self.ylim = None
        self._font_title = font(25 * scale, bold=True)
        self._font_label = font(17 * scale)
        self._font_tick = font(13 * scale)
        self._font_legend = font(13 * scale)

    # ------------------------------------------------------------------ data
    def set_limits(self, xmin, xmax, ymin, ymax, pad=0.02):
        if xmax <= xmin:
            xmax = xmin + 1.0
        if ymax <= ymin:
            ymax = ymin + 1.0
        dx = (xmax - xmin) * pad
        dy = (ymax - ymin) * pad
        self.xlim = (xmin - dx, xmax + dx)
        self.ylim = (ymin - dy, ymax + dy)
        return self

    def _x(self, x):
        x0, x1 = self.xlim
        left, _, right, _ = self.rect
        return left + (x - x0) / (x1 - x0) * (right - left)

    def _y(self, y):
        y0, y1 = self.ylim
        _, top, _, bottom = self.rect
        return bottom - (y - y0) / (y1 - y0) * (bottom - top)

    def _clip_line(self, pts):
        # simple Cohen-Sutherland clipping to the plot rectangle
        left, top, right, bottom = self.rect
        out = []
        for x, y in pts:
            out.append((min(max(x, left), right), min(max(y, top), bottom)))
        return out

    # ------------------------------------------------------------------ axes
    def axes(self, xlabel="", ylabel="", title="", subtitle="",
             xticks=None, yticks=None, xticklabels=None, yticklabels=None,
             grid=True, xfmt=None, yfmt=None):
        if self.xlim is None:
            raise RuntimeError("set_limits() before axes()")
        left, top, right, bottom = self.rect
        if xticks is None:
            xticks = _nice_ticks(self.xlim[0], self.xlim[1])
        if yticks is None:
            yticks = _nice_ticks(self.ylim[0], self.ylim[1])
        if grid:
            for t in xticks:
                x = self._x(t)
                self.draw.line([(x, top), (x, bottom)], fill=(232, 234, 237), width=max(1, int(1 * self.s)))
            for t in yticks:
                y = self._y(t)
                self.draw.line([(left, y), (right, y)], fill=(232, 234, 237), width=max(1, int(1 * self.s)))
        self.draw.rectangle([int(left), int(top), int(right), int(bottom)], outline=(45, 52, 64), width=max(1, int(1.4 * self.s)))

        for i, t in enumerate(xticks):
            x = self._x(t)
            self.draw.line([(x, bottom), (x, bottom + 6 * self.s)], fill=(45, 52, 64),
                           width=max(1, int(1.2 * self.s)))
            label = xticklabels[i] if xticklabels else (xfmt(t) if xfmt else _fmt(t))
            self.draw.text((x, bottom + 10 * self.s), label, font=self._font_tick,
                           fill=(45, 52, 64), anchor="ma")
        for i, t in enumerate(yticks):
            y = self._y(t)
            self.draw.line([(left - 6 * self.s, y), (left, y)], fill=(45, 52, 64),
                           width=max(1, int(1.2 * self.s)))
            label = yticklabels[i] if yticklabels else (yfmt(t) if yfmt else _fmt(t))
            self.draw.text((left - 10 * self.s, y), label, font=self._font_tick,
                           fill=(45, 52, 64), anchor="rm")

        if xlabel:
            self.draw.text(((left + right) / 2, self.H - 24 * self.s), xlabel,
                           font=self._font_label, fill=(30, 35, 45), anchor="mm")
        if ylabel:
            self._rotated_text((26 * self.s, (top + bottom) / 2), ylabel, self._font_label)
        if title:
            self.draw.text(((left + right) / 2, 22 * self.s), title,
                           font=self._font_title, fill=(20, 25, 35), anchor="ma")
            if subtitle:
                self.draw.text(((left + right) / 2, 52 * self.s), subtitle,
                               font=self._font_tick, fill=(90, 98, 110), anchor="ma")
        return self

    def _rotated_text(self, center, text, fnt):
        bbox = self.draw.textbbox((0, 0), text, font=fnt)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        tmp = Image.new("RGBA", (tw + 8, th + 8), (0, 0, 0, 0))
        ImageDraw.Draw(tmp).text((4, 4), text, font=fnt, fill=(30, 35, 45, 255))
        tmp = tmp.rotate(90, expand=True)
        self.img.paste(tmp, (int(center[0] - tmp.width / 2), int(center[1] - tmp.height / 2)), tmp)

    # ------------------------------------------------------------------- draw
    def line(self, xs, ys, color=(31, 58, 147), width=1.6):
        pts = [(self._x(x), self._y(y)) for x, y in zip(xs, ys)]
        if len(pts) >= 2:
            self.draw.line(self._clip_line(pts), fill=color, width=max(1, int(width * self.s)),
                           joint="curve")
        return self

    def gradient_line(self, xs, ys, c0=(31, 58, 147), c1=(192, 57, 43), width=1.2):
        n = min(len(xs), len(ys))
        for i in range(n - 1):
            f = i / max(1, n - 2)
            color = tuple(int(c0[k] + (c1[k] - c0[k]) * f) for k in range(3))
            pts = [(self._x(xs[i]), self._y(ys[i])), (self._x(xs[i + 1]), self._y(ys[i + 1]))]
            self.draw.line(self._clip_line(pts), fill=color, width=max(1, int(width * self.s)))
        return self

    def scatter(self, xs, ys, color=(142, 68, 173), radius=0.7):
        r = max(0.5, radius) * self.s
        for x, y in zip(xs, ys):
            px, py = self._x(x), self._y(y)
            if r <= 0.6 * self.s:
                self.draw.point((px, py), fill=color)
            else:
                self.draw.ellipse([px - r, py - r, px + r, py + r], fill=color)
        return self

    def horizontal(self, y, color=(120, 120, 120), width=1.0, dash=None):
        left, top, right, bottom = self.rect
        py = self._y(y)
        if dash:
            step = dash * self.s
            x = left
            while x < right:
                self.draw.line([(x, py), (min(x + step * 0.55, right), py)], fill=color,
                               width=max(1, int(width * self.s)))
                x += step
        else:
            self.draw.line([(left, py), (right, py)], fill=color, width=max(1, int(width * self.s)))
        return self

    def vertical(self, x, color=(120, 120, 120), width=1.0, dash=None):
        left, top, right, bottom = self.rect
        px = self._x(x)
        if dash:
            step = dash * self.s
            y = top
            while y < bottom:
                self.draw.line([(px, y), (px, min(y + step * 0.55, bottom))], fill=color,
                               width=max(1, int(width * self.s)))
                y += step
        else:
            self.draw.line([(px, top), (px, bottom)], fill=color, width=max(1, int(width * self.s)))
        return self

    def text(self, x, y, message, color=(20, 25, 35), size=14, anchor="la"):
        self.draw.text((self._x(x), self._y(y)), message, font=font(size * self.s),
                       fill=color, anchor=anchor)
        return self

    def legend(self, entries, corner="upper right"):
        # entries: (label, color, style) with style in {"line", "dot"}
        pad = 8 * self.s
        line_h = 21 * self.s
        widths = []
        for label, _, _ in entries:
            bbox = self.draw.textbbox((0, 0), label, font=self._font_legend)
            widths.append(bbox[2] - bbox[0])
        box_w = max(widths) + 48 * self.s
        box_h = line_h * len(entries) + 2 * pad
        left, top, right, bottom = self.rect
        x1 = right - 10 * self.s
        y1 = top + 10 * self.s
        x0 = x1 - box_w
        y0 = y1
        self.draw.rounded_rectangle([x0, y0, x1, y0 + box_h], radius=int(4 * self.s),
                                    fill=(255, 255, 255), outline=(210, 214, 220),
                                    width=max(1, int(1 * self.s)))
        for i, (label, color, style) in enumerate(entries):
            cy = y0 + pad + line_h * i + line_h // 2
            if style == "dot":
                self.draw.ellipse([x0 + 10 * self.s, cy - 4 * self.s,
                                   x0 + 18 * self.s, cy + 4 * self.s], fill=color)
            else:
                self.draw.line([(x0 + 8 * self.s, cy), (x0 + 26 * self.s, cy)], fill=color,
                               width=max(1, int(2 * self.s)))
            self.draw.text((x0 + 32 * self.s, cy), label, font=self._font_legend,
                           fill=(30, 35, 45), anchor="lm")
        return self

    def imshow(self, matrix, cmap=None, vmin=0, vmax=1):
        left, top, right, bottom = self.rect
        h = len(matrix)
        w = len(matrix[0]) if h else 0
        if not cmap:
            cmap = [(255, 255, 255), (31, 58, 147)]
        img = Image.new("RGB", (w, h))
        pix = img.load()
        for j in range(h):
            for i in range(w):
                v = matrix[j][i]
                f = 0.0 if vmax == vmin else (v - vmin) / (vmax - vmin)
                f = min(1.0, max(0.0, f))
                idx = int(f * (len(cmap) - 1))
                pix[i, j] = cmap[idx]
        img = img.resize((right - left, bottom - top), Image.NEAREST)
        self.img.paste(img, (left, top))
        return self

    def save(self, path):
        out = self.img.resize((self.w, self.h), Image.LANCZOS)
        out.save(path)
        return path


def vstack(paths, out, gap=18, background="white"):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    width = max(im.width for im in imgs)
    height = sum(im.height for im in imgs) + gap * (len(imgs) - 1)
    canvas = Image.new("RGB", (width, height), background)
    y = 0
    for im in imgs:
        canvas.paste(im, ((width - im.width) // 2, y))
        y += im.height + gap
    canvas.save(out)
    return out


def hstack(paths, out, gap=18, background="white"):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    width = sum(im.width for im in imgs) + gap * (len(imgs) - 1)
    height = max(im.height for im in imgs)
    canvas = Image.new("RGB", (width, height), background)
    x = 0
    for im in imgs:
        canvas.paste(im, (x, (height - im.height) // 2))
        x += im.width + gap
    canvas.save(out)
    return out
