"""Create small, frame-preserving MP4 previews from the numerical GIF animations."""
import os,sys,argparse,json,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
MEDIA=Path('E:/comsol_runtime/media_tools')
TEMP=Path('E:/comsol_runtime/media_temp')
TEMP.mkdir(parents=True,exist_ok=True)
os.environ['TEMP']=str(TEMP);os.environ['TMP']=str(TEMP)
sys.path.insert(0,str(MEDIA))
import imageio_ffmpeg
from PIL import Image


def make_previews(root):
    executable=imageio_ffmpeg.get_ffmpeg_exe();rows=[]
    for source in sorted(root.glob('**/*_displacement.gif')):
        target=source.with_suffix('.mp4')
        if not target.exists() or target.stat().st_mtime<source.stat().st_mtime:
            print('ENCODE',source,flush=True)
            subprocess.run([executable,'-hide_banner','-loglevel','error','-y','-ignore_loop','1',
                '-i',str(source),'-an','-c:v','libx264','-threads','1','-crf','18','-preset','medium',
                '-pix_fmt','yuv420p','-fps_mode','passthrough','-movflags','+faststart',str(target)],check=True)
        with Image.open(source) as gif:expected=gif.n_frames
        reader=imageio_ffmpeg.read_frames(str(target),output_params=['-vsync','0'])
        metadata=next(reader);actual=sum(1 for _ in reader)
        if actual!=expected:raise RuntimeError(f'Frame count changed: {source}: {expected} -> {actual}')
        rows.append(dict(source=str(source.resolve()),video=str(target.resolve()),frames=actual,
            duration_s=metadata['duration'],source_bytes=source.stat().st_size,video_bytes=target.stat().st_size,
            full_decode_verified=True,numerical_data_changed=False))
        print('VERIFIED',target.name,actual,'frames',target.stat().st_size,'bytes',flush=True)
    (root/'video_manifest.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
    return rows


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8');sys.stderr.reconfigure(encoding='utf-8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'output/transient_matched')
    make_previews(parser.parse_args().output.resolve())
